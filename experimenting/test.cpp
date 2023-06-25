// Vienna Nucleotides:
// 0    1    2    3    4
// N    A    C    G    U

// ContraFold Nucleotides:
// 0    1    2    3    4
// A    C    G    U    N

// Pairings:
// 0     1     2     3     4     5     6     7
// NP    CG    GC    GU    UG    AU    UA    NN

#include <iostream>
#include <string>

using namespace std;

std::string extractNucs(const std::string &input, const std::string &prefix, char delimiter = '_')
{
    std::string out = input.substr(prefix.length());
    // split the trimmed input by delimiter
    int pos = 0;
    for (char c : out)
    {
        if (c == delimiter)
        {
            out.erase(out.begin() + pos);
        }
        pos++;
    }
    return out;
}

int main()
{
    cout << extractVersions("pref_AVA_E_D_DW", "pref") << endl;
    return 0;
}