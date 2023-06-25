#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <map>
#include <chrono>

#include "array_utils.h"

#define VIE_INF 10000000
#define NUCS_NUM 5
#define MAXLOOPSIZE 30
#define NBPAIRS 7

#define NUM_TO_NUC(x) (x == -1 ? -1 : (((x == 4 or x == 5) ? 0 : (x + 1))))
#define NUM_TO_PAIR(x, y) (x == 0 ? (y == 3 ? 5 : 7) : (x == 1 ? (y == 2 ? 1 : 7) : (x == 2 ? (y == 1 ? 2 : (y == 3 ? 3 : 7)) : (x == 3 ? (y == 2 ? 4 : (y == 0 ? 6 : 7)) : 7)))) // 7 is _? or ?_
#define NUC_TO_PAIR(x, y) (x == 1 ? (y == 4 ? 5 : 0) : (x == 2 ? (y == 3 ? 1 : 0) : (x == 3 ? (y == 2 ? 2 : (y == 4 ? 3 : 0)) : (x == 4 ? (y == 3 ? 4 : (y == 1 ? 6 : 0)) : 0))))

double lxc37 = 107.856;
std::string V_NUCS = "ACGUN";
std::string PAIRS[] = {"NP", "CG", "GC", "GU", "UG", "AU", "UA", "NN"};

// prefixes
std::string stackPrefix = "stack_";
std::string hairpinLengthPrefix = "hairpin_length_";
std::string internalLengthPrefix = "internal_length_";
std::string bulgeLengthPrefix = "bulge_length_";
std::string hairpinMismatchPrefix = "terminal_mismatch_hairpin_";
std::string multiMismatchPrefix = "terminal_mismatch_multi_";
std::string externalMismatchPrefix = "terminal_mismatch_external_";
std::string internalMismatchPrefix = "terminal_mismatch_internal_";
std::string internalExplicitPrefix = "internal_explicit_";
std::string dangle5Prefix = "dangle_5_";
std::string dangle3Prefix = "dangle_3_";

// Utility functions
bool decodeEnergyString(std::string &input, const std::string &prefix, char delimiter = '_')
{
    if (input.find(prefix) == std::string::npos)
        return false;

    input = input.substr(prefix.length());
    int pos = 0;
    for (char c : input)
    {
        if (c == delimiter)
            input.erase(input.begin() + pos);
        pos++;
    }
    return true;
}

// Energy Model Class
class EnergyModel
{

private:
    int ML_intern37;
    int ML_closing37;
    int ML_BASE37;
    int MAX_NINIO;
    int ninio37;
    int TerminalAU37; // Outermost pair is AU or GU; also used in tetra_loop triloop
    int **stack37;
    int *hairpin37;
    int *bulge37;
    int *internal_loop37;
    int ***mismatchH37;   // Terminal mismatch for hairpin loop
    int ***mismatchM37;   // Terminal mismatch for multi loop
    int ***mismatchExt37; // Terminal mismatch for external loop
    int ***mismatchI37;   // Terminal mismatch for internal loop
    int ***mismatch1nI37;
    int ***mismatch23I37;
    int ****int11_37;
    int *****int21_37;
    int ******int22_37;
    int **dangle5_37;
    int **dangle3_37;

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

    void updateNucsIdx(std::string &input, int nucsIdxs[])
    {
        for (int i = 0; i < (int)input.size(); i++)
            nucsIdxs[i] = V_NUCS.find(input[i]);
    }

public:
    // takes a file path and loads the energy model
    EnergyModel(const std::string &filepath)
    {

        auto start = std::chrono::high_resolution_clock::now();
        // initialize the arrays
        stack37 = create2DArray(NBPAIRS + 1, NBPAIRS + 1);

        hairpin37 = new int[MAXLOOPSIZE + 1];
        bulge37 = new int[MAXLOOPSIZE + 1];
        internal_loop37 = new int[MAXLOOPSIZE + 1];

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

        int *nucsIdxs = new int[6];
        while (std::getline(file, line))
        {
            std::istringstream is(line);
            std::string key;
            int value;
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
                else
                    std::cout << "Warning: Unknown key in Energy Data: " << key << " -> " << value << std::endl;
            }
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
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
        // #ifdef SPECIAL_HP
        //         if (size == 4 && tetra_hex_tri_index > -1)
        //             return Tetraloop37[tetra_hex_tri_index];
        //         else if (size == 6 && tetra_hex_tri_index > -1)
        //             return Hexaloop37[tetra_hex_tri_index];
        //         else if (size == 3)
        //         {
        //             if (tetra_hex_tri_index > -1)
        //                 return Triloop37[tetra_hex_tri_index];
        //             return (energy + (type > 2 ? TerminalAU37 : 0));
        //         }
        // #endif

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
