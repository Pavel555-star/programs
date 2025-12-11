#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <string>
using namespace std;

string Create_input_from_coordinates(vector<string> species, vector<string> x, vector<string> y, vector<string> z)
    {
    unsigned int i, j;
    unsigned int count_species;
    unsigned int count_layers;
    vector<unsigned int> counts_layers;
    string input = "";
    string ng_1 = "He_";
    string ng_2 = "Ne_";
    string ng_3 = "Ar_";
    string ng_4 = "Kr_";
    string ng_5 = "Xe_";
    string ng_6 = "Rn_";
    string specie;
    size_t found;
    
    string x_string, y_string, z_string;
    
    
    vector<string> elements = {"H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
    
    count_species = species.size();
    counts_layers.reserve(count_species);
    // Check input vectors
    if (count_species != x.size() or count_species != y.size() or count_species != z.size())
        {
        input = "wrong size of coordinate vector";
        return(input);
        }
    // Generate string of atoms
    for (i = 0; i < count_species; i++)
        {
        specie = species[i];
        for (j = 0; j < 118; j++)
            {
            found = specie.find(elements[j]);
            if (found != string::npos)
                {
                input.append(specie + "_");
                count_layers = 1;
                if (j > 2) {
                    input.append(ng_1);
                    count_layers++;
                    }
                if (j > 10) {
                    input.append(ng_2);
                    count_layers++;
                    }
                if (j > 18) {
                    input.append(ng_3);
                    count_layers++;
                    }
                if (j > 36) {
                    input.append(ng_4);
                    count_layers++;
                    }
                if (j > 54) {
                    input.append(ng_5);
                    count_layers++;
                    }
                if (j > 86) {
                    input.append(ng_6);
                    count_layers++;
                    }
                counts_layers.push_back(count_layers);
                break;
                }
            }
        }
    // Append part of coordinates
    for (i = 0; i < count_species; i++)
        {
        x_string = x[i];
        y_string = y[i];
        z_string = z[i];
        count_layers = counts_layers[i];
        for (j = 0; j < count_layers; j++)
            {
            if (i == 0 and j == 0)
                continue;
            
            input.append("[");
            input.append(x_string);
            input.append("_");
            input.append(y_string);
            input.append("_");
            input.append(z_string);
            input.append("]");
            }
        }
    return(input);
    }
    
int main(int argc, char *argv[])
{
vector<string> species ={"Be", "H", "H"};
vector<string> x = {"0.00", "2.60", "0.00"};
vector<string> y = {"0.00", "0.00", "0.00"};
vector<string> z = {"0.00", "0.00", "-2.60"};

string input;
input = Create_input_from_coordinates(species, x, y, z);
cout << input << endl;

return(0);
}
