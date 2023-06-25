#include <iostream>
#include <string>
#include "energy_parameter.h"
#include "intl11.h"
#include "intl21.h"
#include "intl22.h"

int main()
{
    std::string pairs[] = {"NP", "CG", "GC", "GU", "UG", "AU", "UA", "NN"}; // pairs as an array of strings
    std::string nucleotides = "ACGUN";

    // Extract values at the top
    std::cout << "ML_intern37"
              << " " << ML_intern37 << std::endl;
    std::cout << "ML_closing37"
              << " " << ML_closing37 << std::endl;
    std::cout << "ML_BASE37"
              << " " << ML_BASE37 << std::endl;
    std::cout << "MAX_NINIO"
              << " " << MAX_NINIO << std::endl;
    std::cout << "ninio37"
              << " " << ninio37 << std::endl;
    std::cout << "TerminalAU37"
              << " " << TerminalAU37 << std::endl;

    // Triloop values
    for (int i = 0; i < 2; i++)
    {
        int val = Triloop37[i];
        if (val == VIE_INF)
            continue;
        std::cout << "triloop_length_" << i << " " << Triloop37[i] << std::endl;
    }

    // Tetraloop values
    for (int i = 0; i < 16; i++)
    {
        int val = Tetraloop37[i];
        if (val == VIE_INF)
            continue;
        std::cout << "tetraloop_length_" << i << " " << Tetraloop37[i] << std::endl;
    }

    // Hexaloop values
    for (int i = 0; i < 4; i++)
    {
        int val = Hexaloop37[i];
        if (val == VIE_INF)
            continue;
        std::cout << "hexaloop_length_" << i << " " << Hexaloop37[i] << std::endl;
    }

    // Hairpin length values
    for (int i = 0; i < 31; i++)
    {
        int val = hairpin37[i];
        // if (val == VIE_INF)
        //     continue;
        std::cout << "hairpin_length_" << i << " " << hairpin37[i] << std::endl;
    }

    // Bulge length values
    for (int i = 0; i < 31; i++)
    {
        int val = bulge37[i];
        if (val == VIE_INF)
            continue;
        std::cout << "bulge_length_" << i << " " << bulge37[i] << std::endl;
    }

    // Internal loop length values
    for (int i = 0; i < 31; i++)
    {
        int val = internal_loop37[i];
        if (val == VIE_INF)
            continue;
        std::cout << "internal_length_" << i << " " << internal_loop37[i] << std::endl;
    }

    // Stack values
    for (int pair1 = 0; pair1 <= NBPAIRS; pair1++)
    {
        for (int pair2 = 0; pair2 <= NBPAIRS; pair2++)
        {
            int val = stack37[pair1][pair2];
            if (val == VIE_INF)
                continue;
            std::cout << "stack_" << pairs[pair1] << pairs[pair2] << " " << stack37[pair1][pair2] << std::endl;
        }
    }

    // Terminal mismatch internal loop values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        std::string pair_str = pairs[pair]; // get the current pair as a string
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                int val = mismatchI37[pair][i][j];
                if (val == VIE_INF)
                    continue;
                std::cout << "terminal_mismatch_internal_" << pair_str << "_" << nucleotides[i] << "_" << nucleotides[j]
                          << " " << mismatchI37[pair][i][j] << std::endl;
            }
        }
    }

    // Terminal mismatch for hairpin loop values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        std::string pair_str = pairs[pair]; // get the current pair as a string
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                // int val = mismatchH37[pair][i][j];
                // if (val == VIE_INF)
                //     continue;
                std::cout << "terminal_mismatch_hairpin_" << pair_str << "_" << nucleotides[i] << "_" << nucleotides[j]
                          << " " << mismatchH37[pair][i][j] << std::endl;
            }
        }
    }

    // Terminal Mismatch for multi Loop values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        std::string pair_str = pairs[pair]; // get the current pair as a string
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                int val = mismatchM37[pair][i][j];
                if (val == VIE_INF)
                    continue;
                std::cout << "terminal_mismatch_multi_" << pair_str << "_" << nucleotides[i] << "_" << nucleotides[j]
                          << " " << mismatchM37[pair][i][j] << std::endl;
            }
        }
    }

    // Terminal Mismatch for external loop values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        std::string pair_str = pairs[pair]; // get the current pair as a string
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                int val = mismatchM37[pair][i][j];
                if (val == VIE_INF)
                    continue;
                std::cout << "terminal_mismatch_external_" << pair_str << "_" << nucleotides[i] << "_" << nucleotides[j]
                          << " " << mismatchExt37[pair][i][j] << std::endl;
            }
        }
    }

    // Terminal Mismatch for internal loop 1-n values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                std::string pair_str = pairs[pair]; // get the current pair as a string
                int val = mismatch1nI37[pair][i][j];
                if (val == VIE_INF)
                    continue;
                std::cout << "terminal_mismatch_internal_1n_" << pair_str << "_" << nucleotides[i] << "_"
                          << nucleotides[j] << " " << mismatch1nI37[pair][i][j] << std::endl;
            }
        }
    }

    // Terminal Mismatch for internal loop 2-3 values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                std::string pair_str = pairs[pair]; // get the current pair as a string
                int val = mismatch23I37[pair][i][j];
                if (val == VIE_INF)
                    continue;
                std::cout << "terminal_mismatch_internal_23_" << pair_str << "_" << nucleotides[i] << "_"
                          << nucleotides[j] << " " << mismatch23I37[pair][i][j] << std::endl;
            }
        }
    }

    // Extract the 1-1 internal loop values
    // 4 Dimensions: pair, pair, nucletoid, nucleotide
    for (int pair1 = 0; pair1 <= NBPAIRS; pair1++)
    {
        for (int pair2 = 0; pair2 <= NBPAIRS; pair2++)
        {
            std::string pair1_str = pairs[pair1]; // get the current pair as a string
            std::string pair2_str = pairs[pair2]; // get the current pair as a string
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    int val = int11_37[pair1][pair2][i][j];
                    if (val == VIE_INF)
                        continue;
                    std::cout << "internal_explicit_1_1_" << pair1_str << "_" << pair2_str << "_" << nucleotides[i] << "_"
                              << nucleotides[j] << " " << int11_37[pair1][pair2][i][j] << std::endl;
                }
            }
        }
    }

    // Extract the 2-1 internal loop values
    // 5 Dimensions: pair, pair, nucletoid, nucleotide, nucleotide
    for (int pair1 = 0; pair1 <= NBPAIRS; pair1++)
    {
        for (int pair2 = 0; pair2 <= NBPAIRS; pair2++)
        {
            std::string pair1_str = pairs[pair1]; // get the current pair as a string
            std::string pair2_str = pairs[pair2]; // get the current pair as a string
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    for (int k = 0; k < 5; k++)
                    {
                        int val = int21_37[pair1][pair2][i][j][k];
                        if (val == VIE_INF)
                            continue;
                        std::cout << "internal_explicit_2_1_" << pair1_str << "_" << pair2_str << "_" << nucleotides[i]
                                  << "_" << nucleotides[j] << "_" << nucleotides[k] << " "
                                  << int21_37[pair1][pair2][i][j][k] << std::endl;
                    }
                }
            }
        }
    }

    // Extract the 2-2 internal loop values
    // 6 Dimensions: pair, pair, nucletoid, nucleotide, nucleotide, nucleotide
    for (int pair1 = 0; pair1 <= NBPAIRS; pair1++)
    {
        for (int pair2 = 0; pair2 <= NBPAIRS; pair2++)
        {
            std::string pair1_str = pairs[pair1]; // get the current pair as a string
            std::string pair2_str = pairs[pair2]; // get the current pair as a string
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    for (int k = 0; k < 5; k++)
                    {
                        for (int l = 0; l < 5; l++)
                        {
                            int val = int22_37[pair1][pair2][i][j][k][l];
                            if (val == VIE_INF)
                                continue;
                            std::cout << "internal_explicit_2_2_" << pair1_str << "_" << pair2_str << "_"
                                      << nucleotides[i] << "_" << nucleotides[j] << "_" << nucleotides[k] << "_"
                                      << nucleotides[l] << " " << int22_37[pair1][pair2][i][j][k][l] << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Dangle5 values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        std::string pair_str = pairs[pair]; // get the current pair as a string
        for (int i = 0; i < 5; i++)
        {
            int val = dangle5_37[pair][i];
            if (val == VIE_INF)
                continue;
            std::cout << "dangle_5_" << pair_str << "_" << nucleotides[i] << " " << dangle5_37[pair][i] << std::endl;
        }
    }

    // Dangle3 values
    for (int pair = 0; pair <= NBPAIRS; pair++)
    {
        std::string pair_str = pairs[pair]; // get the current pair as a string
        for (int i = 0; i < 5; i++)
        {
            int val = dangle3_37[pair][i];
            if (val == VIE_INF)
                continue;
            std::cout << "dangle_3_" << pair_str << "_" << nucleotides[i] << " " << dangle3_37[pair][i] << std::endl;
        }
    }

    return 0;
}
