#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

#include "energy_model.h"
#include "array_utils.h"

// Include OpenMP for parallelization (for parsing energy text file) if available
#ifdef _OPENMP
#include <omp.h>
#endif

// Nucleotides
const std::string C_NUCS = "ACGUN"; // Contrafold
const std::string V_NUCS = "NACGU"; // Vienna
const std::string PAIRS[] = {"NP", "CG", "GC", "GU", "UG", "AU", "UA", "NN"};

// Parsing prefixes
const std::string stackPrefix = "stack_";
const std::string triloopPrefix = "triloop_length_";
const std::string tetraloopPrefix = "tetraloop_length_";
const std::string hexaloopPrefix = "hexaloop_length_";
const std::string hairpinLengthPrefix = "hairpin_length_";
const std::string internalLengthPrefix = "internal_length_";
const std::string bulgeLengthPrefix = "bulge_length_";
const std::string hairpinMismatchPrefix = "terminal_mismatch_hairpin_";
const std::string multiMismatchPrefix = "terminal_mismatch_multi_";
const std::string externalMismatchPrefix = "terminal_mismatch_external_";
const std::string internalMismatchPrefix = "terminal_mismatch_internal_";
const std::string internalExplicitPrefix = "internal_explicit_";
const std::string dangle5Prefix = "dangle_5_";
const std::string dangle3Prefix = "dangle_3_";

// Energy parameters
int ML_intern37;      // ???
int ML_closing37;     // ???
int ML_BASE37;        // ???
int MAX_NINIO;        // ???
int ninio37;          // ???
int TerminalAU37;     // Outermost pair is AU or GU; also used in tetra_loop triloop
int *Triloop37;       // Triloop energies
int *Tetraloop37;     // Tetraloop energies
int *Hexaloop37;      // Hexaloop energies
int **stack37;        // Stacking energies
int *hairpin37;       // Hairpin loop energies (based on length)
int *bulge37;         // Bulge loop energies (based on length)
int *internal_loop37; // Internal loop energies (based on length)
int ***mismatchH37;   // Terminal mismatch energies for hairpin loop
int ***mismatchM37;   // Terminal mismatch energies for multi loop
int ***mismatchExt37; // Terminal mismatch energies for external loop
int ***mismatchI37;   // Terminal mismatch energies for internal loop
int ***mismatch1nI37; // Terminal mismatch energies for internal (1 x N) loop
int ***mismatch23I37; // Terminal mismatch energies for internal (2 x 3) loop
int ****int11_37;     // Terminal mismatch energies for internal (1 x 1) loop
int *****int21_37;    // Terminal mismatch energies for internal (2 x 1) loop
int ******int22_37;   // Terminal mismatch energies for internal (2 x 2) loop
int **dangle5_37;     // Dangle energies for 5' end
int **dangle3_37;     // Dangle energies for 3' end

// Private helper function used for parsing energy data
bool decodeEnergyString(std::string &input, const std::string &prefix, char delimiter = '_')
{
    if (input.find(prefix) == std::string::npos)
        return false;

    input = input.substr(prefix.length());
    input.erase(std::remove(input.begin(), input.end(), delimiter), input.end());
    return true;
}

// Private helper function used for parsing energy data
void updateNucsIdx(std::string &input, int nucsIdxs[])
{
    for (int i = 0; i < (int)input.size(); i++)
        nucsIdxs[i] = C_NUCS.find(input[i]);
}

// Parse energy data from text file, and store in global variables. Publically accessible.
void parseEnergyData(const std::string &filepath, bool verbose)
{
    auto start = std::chrono::high_resolution_clock::now();

    // Initialize the arrays
    Triloop37 = new int[2];
    Tetraloop37 = new int[16];
    Hexaloop37 = new int[4];

    hairpin37 = new int[MAXLOOPSIZE + 1];
    bulge37 = new int[MAXLOOPSIZE + 1];
    internal_loop37 = new int[MAXLOOPSIZE + 1];

    stack37 = create2DArray(NBPAIRS + 1, NBPAIRS + 1);
    mismatchH37 = create3DArray(NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    mismatchM37 = create3DArray(NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    mismatchExt37 = create3DArray(NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    mismatchI37 = create3DArray(NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    mismatch1nI37 = create3DArray(NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    mismatch23I37 = create3DArray(NBPAIRS + 1, NUCS_NUM, NUCS_NUM);

    int11_37 = create4DArray(NBPAIRS + 1, NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    int21_37 = create5DArray(NBPAIRS + 1, NBPAIRS + 1, NUCS_NUM, NUCS_NUM, NUCS_NUM);
    int22_37 = create6DArray(NBPAIRS + 1, NBPAIRS + 1, NUCS_NUM, NUCS_NUM, NUCS_NUM, NUCS_NUM);

    dangle5_37 = create2DArray(NBPAIRS + 1, NUCS_NUM);
    dangle3_37 = create2DArray(NBPAIRS + 1, NUCS_NUM);

    std::ifstream file(filepath);
    std::string line;

    std::vector<std::string> lines;

    // Serial reading of the file
    while (std::getline(file, line))
        lines.push_back(line);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < lines.size(); i++)
    {
        int *nucsIdxs = new int[10], value;
        std::istringstream is(lines[i]);
        std::string key;
        if (is >> key >> value)
        {
            if (decodeEnergyString(key, "ML_intern"))
                ML_intern37 = value;

            else if (decodeEnergyString(key, "ML_closing"))
                ML_closing37 = value;

            else if (decodeEnergyString(key, "ML_BASE"))
                ML_BASE37 = value;

            else if (decodeEnergyString(key, "MAX_NINIO"))
                MAX_NINIO = value;

            else if (decodeEnergyString(key, "ninio"))
                ninio37 = value;

            else if (decodeEnergyString(key, "TerminalAU"))
                TerminalAU37 = value;

            else if (decodeEnergyString(key, triloopPrefix))
                Triloop37[std::stoi(key)] = value;

            else if (decodeEnergyString(key, tetraloopPrefix))
                Tetraloop37[std::stoi(key)] = value;

            else if (decodeEnergyString(key, hexaloopPrefix))
                Hexaloop37[std::stoi(key)] = value;

            else if (decodeEnergyString(key, hairpinLengthPrefix))
                hairpin37[std::stoi(key)] = value;

            else if (decodeEnergyString(key, internalLengthPrefix))
                internal_loop37[std::stoi(key)] = value;

            else if (decodeEnergyString(key, bulgeLengthPrefix))
                bulge37[std::stoi(key)] = value;

            else if (decodeEnergyString(key, stackPrefix))
            {
                updateNucsIdx(key, nucsIdxs);
                stack37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][NUM_TO_PAIR(nucsIdxs[2], nucsIdxs[3])] = value;
            }

            else if (decodeEnergyString(key, hairpinMismatchPrefix))
            {
                updateNucsIdx(key, nucsIdxs);
                mismatchH37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][nucsIdxs[2]][nucsIdxs[3]] = value;
            }

            else if (decodeEnergyString(key, multiMismatchPrefix))
            {
                updateNucsIdx(key, nucsIdxs);
                mismatchM37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][nucsIdxs[2]][nucsIdxs[3]] = value;
            }

            else if (decodeEnergyString(key, externalMismatchPrefix))
            {
                updateNucsIdx(key, nucsIdxs);
                mismatchExt37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][nucsIdxs[2]][nucsIdxs[3]] = value;
            }

            else if (decodeEnergyString(key, dangle5Prefix))
            {
                updateNucsIdx(key, nucsIdxs);
                dangle5_37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][nucsIdxs[2]] = value;
            }

            else if (decodeEnergyString(key, dangle3Prefix))
            {
                updateNucsIdx(key, nucsIdxs);
                dangle3_37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][nucsIdxs[2]] = value;
            }

            else if (decodeEnergyString(key, internalMismatchPrefix))
            {
                updateNucsIdx(key, nucsIdxs);
                if (key[0] == '1') // Corresponds to mismatch1nI37
                    mismatch1nI37[NUM_TO_PAIR(nucsIdxs[2], nucsIdxs[3])][nucsIdxs[4]][nucsIdxs[5]] = value;
                else if (key[0] == '2') // Corresponds to mismatch23I37
                    mismatch23I37[NUM_TO_PAIR(nucsIdxs[2], nucsIdxs[3])][nucsIdxs[4]][nucsIdxs[5]] = value;
                else // Corresponds to mismatchI37
                    mismatchI37[NUM_TO_PAIR(nucsIdxs[0], nucsIdxs[1])][nucsIdxs[2]][nucsIdxs[3]] = value;
            }

            else if (decodeEnergyString(key, internalExplicitPrefix))
            {
                updateNucsIdx(key, nucsIdxs);
                if (key[0] == '1') // Corresponds to int11_37
                    int11_37[NUM_TO_PAIR(nucsIdxs[2], nucsIdxs[3])][NUM_TO_PAIR(nucsIdxs[4], nucsIdxs[5])][nucsIdxs[6]][nucsIdxs[7]] = value;
                else
                {
                    if (key[1] == '1') // Corresponds to int21_37
                        int21_37[NUM_TO_PAIR(nucsIdxs[2], nucsIdxs[3])][NUM_TO_PAIR(nucsIdxs[4], nucsIdxs[5])][nucsIdxs[6]][nucsIdxs[7]][nucsIdxs[8]] = value;
                    else // Corresponds to int22_37
                        int22_37[NUM_TO_PAIR(nucsIdxs[2], nucsIdxs[3])][NUM_TO_PAIR(nucsIdxs[4], nucsIdxs[5])][nucsIdxs[6]][nucsIdxs[7]][nucsIdxs[8]][nucsIdxs[9]] = value;
                }
            }
            else if (verbose)
            {
                std::cout << "Warning: Unknown key in Energy Data: " << key << " -> " << value << std::endl;
            }

            // Free memory
            delete[] nucsIdxs;
        }
    }

    lines.clear();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    if (verbose)
        std::cout << "Time to load energy model: " << elapsed.count() << " s\n";
}

// To free memory and avoid memory leaks. Publically accessible.
void freeMemory()
{
    delete2DArray(stack37, NBPAIRS + 1);
    delete[] hairpin37;
    delete[] bulge37;
    delete[] internal_loop37;
    delete3DArray(mismatchH37, NBPAIRS + 1, NUCS_NUM);
    delete3DArray(mismatchM37, NBPAIRS + 1, NUCS_NUM);
    delete3DArray(mismatchExt37, NBPAIRS + 1, NUCS_NUM);
    delete3DArray(mismatchI37, NBPAIRS + 1, NUCS_NUM);
    delete3DArray(mismatch1nI37, NBPAIRS + 1, NUCS_NUM);
    delete3DArray(mismatch23I37, NBPAIRS + 1, NUCS_NUM);
    delete4DArray(int11_37, NBPAIRS + 1, NBPAIRS + 1, NUCS_NUM);
    delete5DArray(int21_37, NBPAIRS + 1, NBPAIRS + 1, NUCS_NUM, NUCS_NUM);
    delete6DArray(int22_37, NBPAIRS + 1, NBPAIRS + 1, NUCS_NUM, NUCS_NUM, NUCS_NUM);
    delete2DArray(dangle5_37, NBPAIRS + 1);
    delete2DArray(dangle3_37, NBPAIRS + 1);
}

////////////////////////////////////////////////// Array Helpers //////////////////////////////////////////////////

int **create2DArray(int d1, int d2)
{
    int **arr = new int *[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = new int[d2];
        std::fill_n(arr[i], d2, VIE_INF);
    }
    return arr;
}

int ***create3DArray(int d1, int d2, int d3)
{
    int ***arr = new int **[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create2DArray(d2, d3);
    }
    return arr;
}

int ****create4DArray(int d1, int d2, int d3, int d4)
{
    int ****arr = new int ***[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create3DArray(d2, d3, d4);
    }
    return arr;
}

int *****create5DArray(int d1, int d2, int d3, int d4, int d5)
{
    int *****arr = new int ****[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create4DArray(d2, d3, d4, d5);
    }
    return arr;
}

int ******create6DArray(int d1, int d2, int d3, int d4, int d5, int d6)
{
    int ******arr = new int *****[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create5DArray(d2, d3, d4, d5, d6);
    }
    return arr;
}

void delete2DArray(int **arr, int d1)
{
    for (int i = 0; i < d1; ++i)
    {
        delete[] arr[i];
    }
    delete[] arr;
}

void delete3DArray(int ***arr, int d1, int d2)
{
    for (int i = 0; i < d1; ++i)
    {
        delete2DArray(arr[i], d2);
    }
    delete[] arr;
}

void delete4DArray(int ****arr, int d1, int d2, int d3)
{
    for (int i = 0; i < d1; ++i)
    {
        delete3DArray(arr[i], d2, d3);
    }
    delete[] arr;
}

void delete5DArray(int *****arr, int d1, int d2, int d3, int d4)
{
    for (int i = 0; i < d1; ++i)
    {
        delete4DArray(arr[i], d2, d3, d4);
    }
    delete[] arr;
}

void delete6DArray(int ******arr, int d1, int d2, int d3, int d4, int d5)
{
    for (int i = 0; i < d1; ++i)
    {
        delete5DArray(arr[i], d2, d3, d4, d5);
    }
    delete[] arr;
}
