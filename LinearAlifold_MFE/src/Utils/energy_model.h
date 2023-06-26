#ifndef ENERGY_MODEL_H
#define ENERGY_MODEL_H

#include <string>

#define VIE_INF 10000000
#define NUCS_NUM 5
#define NBPAIRS 7
#define MAXLOOPSIZE 30
#define lxc37 107.856

#define NUM_TO_NUC(x) (x == -1 ? -1 : (((x == 4 or x == 5) ? 0 : (x + 1))))
#define NUM_TO_PAIR(x, y) (x == 0 ? (y == 3 ? 5 : 7) : (x == 1 ? (y == 2 ? 1 : 7) : (x == 2 ? (y == 1 ? 2 : (y == 3 ? 3 : 7)) : (x == 3 ? (y == 2 ? 4 : (y == 0 ? 6 : 7)) : 7))))
#define NUC_TO_PAIR(x, y) (x == 1 ? (y == 4 ? 5 : 0) : (x == 2 ? (y == 3 ? 1 : 0) : (x == 3 ? (y == 2 ? 2 : (y == 4 ? 3 : 0)) : (x == 4 ? (y == 3 ? 4 : (y == 1 ? 6 : 0)) : 0))))

// Variables for storing energy data.
extern int ML_intern37;      // ???
extern int ML_closing37;     // ???
extern int ML_BASE37;        // ???
extern int MAX_NINIO;        // ???
extern int ninio37;          // ???
extern int TerminalAU37;     // Outermost pair is AU or GU; also used in tetra_loop triloop
extern int *Triloop37;       // Triloop energies
extern int *Tetraloop37;     // Tetraloop energies
extern int *Hexaloop37;      // Hexaloop energies
extern int **stack37;        // Stacking energies
extern int *hairpin37;       // Hairpin loop energies (based on length)
extern int *bulge37;         // Bulge loop energies (based on length)
extern int *internal_loop37; // Internal loop energies (based on length)
extern int ***mismatchH37;   // Terminal mismatch energies for hairpin loop
extern int ***mismatchM37;   // Terminal mismatch energies for multi loop
extern int ***mismatchExt37; // Terminal mismatch energies for external loop
extern int ***mismatchI37;   // Terminal mismatch energies for internal loop
extern int ***mismatch1nI37; // Terminal mismatch energies for internal (1 x N) loop
extern int ***mismatch23I37; // Terminal mismatch energies for internal (2 x 3) loop
extern int ****int11_37;     // Terminal mismatch energies for internal (1 x 1) loop
extern int *****int21_37;    // Terminal mismatch energies for internal (2 x 1) loop
extern int ******int22_37;   // Terminal mismatch energies for internal (2 x 2) loop
extern int **dangle5_37;     // Dangle energies for 5' end
extern int **dangle3_37;     // Dangle energies for 3' end

// Functions for parsing energy data and freeing memory.
void parseEnergyData(const std::string &filepath, bool verbose = true);
void freeMemory();

// Energy model functions. Publicly accessible.
inline int score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int tetra_hex_tri_index = -1)
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

inline int score_single_alifold(int n1, int n2, int type, int type_2,
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
inline int E_MLstem(int type, int si1, int sj1)
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

inline int score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len)
{
    int tt = NUM_TO_PAIR(nucj, nuci);
    int si1 = NUM_TO_NUC(nuci1);
    int sj1 = NUM_TO_NUC(nucj_1);

    return E_MLstem(tt, sj1, si1) + ML_closing37;
}

inline int score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len)
{
    // int p = i;
    // int q = k;
    int tt = NUM_TO_PAIR(nuci, nuck);
    int sp1 = NUM_TO_NUC(nuci_1);
    int sq1 = NUM_TO_NUC(nuck1);

    return E_MLstem(tt, sp1, sq1);
}

// exterior_loop
inline int score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len)
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

inline int score_multi_unpaired(int i, int j)
{
    return 0;
}

inline int score_external_unpaired(int i, int j)
{
    return 0;
}

#endif // ENERGY_MODEL_H