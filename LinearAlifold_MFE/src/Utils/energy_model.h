#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include "array_utils.h"

// Include OpenMP for parallelization (for parsing energy text file) if available
#ifdef _OPENMP
#include <omp.h>
#endif

#define VIE_INF 10000000
#define NUCS_NUM 5
#define MAXLOOPSIZE 30
#define NBPAIRS 7
#define lxc37 107.856

#define NUM_TO_NUC(x) (x == -1 ? -1 : (((x == 4 or x == 5) ? 0 : (x + 1))))
#define NUM_TO_PAIR(x, y) (x == 0 ? (y == 3 ? 5 : 7) : (x == 1 ? (y == 2 ? 1 : 7) : (x == 2 ? (y == 1 ? 2 : (y == 3 ? 3 : 7)) : (x == 3 ? (y == 2 ? 4 : (y == 0 ? 6 : 7)) : 7)))) // 7 is _? or ?_
#define NUC_TO_PAIR(x, y) (x == 1 ? (y == 4 ? 5 : 0) : (x == 2 ? (y == 3 ? 1 : 0) : (x == 3 ? (y == 2 ? 2 : (y == 4 ? 3 : 0)) : (x == 4 ? (y == 3 ? 4 : (y == 1 ? 6 : 0)) : 0))))

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

// Energy Model Class
class EnergyModel
{

private:
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

    bool decodeEnergyString(std::string &input, const std::string &prefix, char delimiter = '_')
    {
        if (input.find(prefix) == std::string::npos)
            return false;

        input = input.substr(prefix.length());
        input.erase(std::remove(input.begin(), input.end(), delimiter), input.end());
        return true;
    }

    void updateNucsIdx(std::string &input, int nucsIdxs[])
    {
        for (int i = 0; i < (int)input.size(); i++)
            nucsIdxs[i] = C_NUCS.find(input[i]);
    }

public:
    // Takes a file path and loads the energy model
    EnergyModel(const std::string &filepath, bool verbose = false)
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

    // Destructor to free memory and avoid memory leaks
    ~EnergyModel()
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

    int score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int tetra_hex_tri_index = -1)
    {
        int size = j - i - 1;
        int type = NUM_TO_PAIR(nuci, nucj);
        int si1 = NUM_TO_NUC(nuci1);
        int sj1 = NUM_TO_NUC(nucj_1);
        int energy;

        if (size <= 30)
            energy = hairpin37[size];
        else
            energy = hairpin37[30] + (int)(lxc37 * log((size) / 30.));
        if (size < 3)
            return energy;
/* should only be the case when folding alignments */
#ifdef SPECIAL_HP
        if (size == 4 && tetra_hex_tri_index > -1)
            return Tetraloop37[tetra_hex_tri_index];
        else if (size == 6 && tetra_hex_tri_index > -1)
            return Hexaloop37[tetra_hex_tri_index];
        else if (size == 3)
        {
            if (tetra_hex_tri_index > -1)
                return Triloop37[tetra_hex_tri_index];
            return (energy + (type > 2 ? TerminalAU37 : 0));
        }
#endif

        energy += mismatchH37[type][si1][sj1];
        return energy;
    };

    int score_single_alifold(int n1, int n2, int type, int type_2,
                             int nuci1, int nucj_1,
                             int nucp_1, int nucq1)
    {

        int si1 = NUM_TO_NUC(nuci1);
        int sj1 = NUM_TO_NUC(nucj_1);
        int sp1 = NUM_TO_NUC(nucp_1);
        int sq1 = NUM_TO_NUC(nucq1);
        int nl, ns, u, energy;
        energy = 0;

        if (n1 > n2)
        {
            nl = n1;
            ns = n2;
        }
        else
        {
            nl = n2;
            ns = n1;
        }

        if (nl == 0)
            return stack37[type][type_2]; /* stack */

        if (ns == 0)
        { /* bulge */
            energy = (nl <= MAXLOOPSIZE) ? bulge37[nl] : (bulge37[30] + (int)(lxc37 * log(nl / 30.)));
            if (nl == 1)
                energy += stack37[type][type_2];
            else
            {
                if (type > 2)
                    energy += TerminalAU37;
                if (type_2 > 2)
                    energy += TerminalAU37;
            }
            return energy;
        }
        else
        { /* interior loop */
            if (ns == 1)
            {
                if (nl == 1) /* 1x1 loop */
                    return int11_37[type][type_2][si1][sj1];
                if (nl == 2)
                { /* 2x1 loop */
                    if (n1 == 1)
                        energy = int21_37[type][type_2][si1][sq1][sj1];
                    else
                        energy = int21_37[type_2][type][sq1][si1][sp1];
                    return energy;
                }
                else
                { /* 1xn loop */
                    energy = (nl + 1 <= MAXLOOPSIZE) ? (internal_loop37[nl + 1]) : (internal_loop37[30] + (int)(lxc37 * log((nl + 1) / 30.)));
                    energy += std::min(MAX_NINIO, (nl - ns) * ninio37);
                    energy += mismatch1nI37[type][si1][sj1] + mismatch1nI37[type_2][sq1][sp1];
                    return energy;
                }
            }
            else if (ns == 2)
            {
                if (nl == 2)
                { /* 2x2 loop */
                    return int22_37[type][type_2][si1][sp1][sq1][sj1];
                }
                else if (nl == 3)
                { /* 2x3 loop */
                    energy = internal_loop37[5] + ninio37;
                    energy += mismatch23I37[type][si1][sj1] + mismatch23I37[type_2][sq1][sp1];
                    return energy;
                }
            }
            { /* generic interior loop (no else here!)*/
                u = nl + ns;
                energy = (u <= MAXLOOPSIZE) ? (internal_loop37[u]) : (internal_loop37[30] + (int)(lxc37 * log((u) / 30.)));

                energy += std::min(MAX_NINIO, (nl - ns) * ninio37);

                energy += mismatchI37[type][si1][sj1] + mismatchI37[type_2][sq1][sp1];
            }
        }
        return energy;
    }

    // multi_loop
    int E_MLstem(int type, int si1, int sj1)
    {
        int energy = 0;

        if (si1 >= 0 && sj1 >= 0)
        {
            energy += mismatchM37[type][si1][sj1];
        }
        else if (si1 >= 0)
        {
            energy += dangle5_37[type][si1];
        }
        else if (sj1 >= 0)
        {
            energy += dangle3_37[type][sj1];
        }

        if (type > 2)
        {
            energy += TerminalAU37;
        }

        energy += ML_intern37;

        return energy;
    }

    int score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len)
    {
        int tt = NUM_TO_PAIR(nucj, nuci);
        int si1 = NUM_TO_NUC(nuci1);
        int sj1 = NUM_TO_NUC(nucj_1);

        return E_MLstem(tt, sj1, si1) + ML_closing37;
    }

    int score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len)
    {
        // int p = i;
        // int q = k;
        int tt = NUM_TO_PAIR(nuci, nuck);
        int sp1 = NUM_TO_NUC(nuci_1);
        int sq1 = NUM_TO_NUC(nuck1);

        return E_MLstem(tt, sp1, sq1);
    }

    // exterior_loop
    int score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len)
    {
        int type = NUM_TO_PAIR(nuci, nucj);
        int si1 = NUM_TO_NUC(nuci_1);
        int sj1 = NUM_TO_NUC(nucj1);
        int energy = 0;

        if (si1 >= 0 && sj1 >= 0)
        {
            energy += mismatchExt37[type][si1][sj1];
        }
        else if (si1 >= 0)
        {
            energy += dangle5_37[type][si1];
        }
        else if (sj1 >= 0)
        {
            energy += dangle3_37[type][sj1];
        }

        if (type > 2)
            energy += TerminalAU37;
        return energy;
    }

    int score_multi_unpaired(int i, int j)
    {
        return 0;
    }

    int score_external_unpaired(int i, int j)
    {
        return 0;
    }
};
