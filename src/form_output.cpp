#include "nuclide_class.h"

namespace openbps {

    std::vector<std::pair<int, std::string>> Chain::form_idx_name() {
        std::vector<std::pair<int, std::string>> out;
        std::pair<int, std::string> elem;

        for (int i = 0; i < nuclides.size(); i++) {
            elem.first = i;
            elem.second = nuclides[i].name;
            out.push_back(elem);
          //  std::cout << "idx" << elem.first << "name" << elem.second << std::endl;
        }
        return out;
    }

    std::vector<std::pair<int, double>> Chain::form_idx_lambda() {
        std::vector<std::pair<int, double>> out;
        std::pair<int, double> elem;

        for (int i = 0; i < nuclides.size(); i++) {
            elem.first = i;
            elem.second = nuclides[i].half_life;
            out.push_back(elem);
            // std::cout << "idx" << elem.first << "Lamda " << elem.second << std::endl;
        }
        return out;
    }

    std::vector<std::vector<double>> Chain::form_idx_decay() {
        std::vector<std::vector<double>> out;
        std::vector<double> elem;

        for (int i = 0; i < nuclides.size(); i++) {
            for (int j = 0; j < nuclides[i].decay_arr.size(); j++) {
                elem.push_back(i);
                elem.push_back(name_idx[nuclides[i].decay_arr[j].target]);
                elem.push_back(nuclides[i].decay_arr[j].branching_ratio);
                out.push_back(elem);
                elem.clear();
            }
        }
        return (out);
    }
/*
 * to test the following parser use:
 *
 *
    auto react_map = chain.form_reaction();
    for(int i = 0; i < react_map["(n,gamma)"].size(); i++)
    {
        std::cout << react_map["(n,gamma)"][i].first << "  " << react_map["(n,gamma)"][i].second << "  "<< 1 << "  " <<std::endl;
    }

    output:

0  1  1
1  2  1
7  8  1
17  18  1
18  19  1
28  29  1
40  41  1
41  42  1

...
 */
    std::map<std::string, std::vector<std::pair<int, int>>> Chain::form_reaction() {

        std::map<std::string, std::vector<std::pair<int, int>>> out;
        std::vector<std::pair<int, int>> n_2n;
        std::vector<std::pair<int, int>> n_3n;
        std::vector<std::pair<int, int>> n_a;
        std::vector<std::pair<int, int>> n_gamma;
        std::vector<std::pair<int, int>> n_p;
        std::pair<int, int> elem;

        out.insert({"(n,2n)", n_2n});
        out.insert({"(n,3n)", n_3n});
        out.insert({"(n,gamma)", n_gamma});
        out.insert({"(n,a)", n_a});
        out.insert({"(n,p)", n_p});

        for (int i = 0; i < nuclides.size(); i++) {
            for (int j = 0; j < nuclides[i].reaction_arr.size(); j++) {
                elem.first = i;
                elem.second = name_idx[nuclides[i].reaction_arr[j].target];
                out[nuclides[i].reaction_arr[j].type].push_back(elem);
            }
        }
        return (out);
    }

/*
 * To test this parser (only Pu239 have 2000000.0):
    auto map = chain.form_yield_map();
    for(int i = 0; i < map["2000000.0"].size(); i++)
    {
        std::cout << map["2000000.0"][i][0] << "  " << map["2000000.0"][i][1] << "  "<< map["2000000.0"][i][2] << "  " <<std::endl;
    }
output
...
3557  1346  0
3557  1350  0
3557  1351  0
3557  1352  6.12916e-14
3557  1353  1.65977e-13
3557  1354  2.47967e-12

and so on

*/

    std::map<std::string, std::vector<std::vector<double>>>
            Chain::form_yield_map(){

        std::map<std::string, std::vector<std::vector<double>>> out;
        std::vector<std::vector<double>> en_1;
        std::vector<std::vector<double>> en_2;
        std::vector<std::vector<double>> en_3;
        std::vector<std::vector<double>> en_4;
        std::vector<double> elem;

        out.insert({"0.0253", en_1});
        out.insert({"500000.0", en_2});
        out.insert({"2000000.0", en_3});
        out.insert({"14000000.0", en_4});

        for(int i = 0; i < nuclides.size(); i++) {

            for(int j = 0 ; j < nuclides[i].nfy.energies.size(); j++) {

                auto n_map = nuclides[i].nfy.yield_arr[j].product_data;
                for(auto& item : n_map)
                {
                    elem.push_back(i);
                    elem.push_back(name_idx[item.first]);
                    elem.push_back(item.second);
                    //std::cout<<elem[0]<<"   "<<elem[1]<<"   " <<elem[2]<<std::endl;
                    if(nuclides[i].nfy.yield_arr[j].energy == 0.0253)
                        out["0.0253"].push_back(elem);
                    if(nuclides[i].nfy.yield_arr[j].energy == 500000.0)
                        out["500000.0"].push_back(elem);
                    if(nuclides[i].nfy.yield_arr[j].energy == 2000000.0)
                        out["2000000.0"].push_back(elem);
                    if(nuclides[i].nfy.yield_arr[j].energy == 14000000.0)
                        out["14000000.0"].push_back(elem);
                    elem.clear();
                }
            }
        }
        return out;
    }

} // namespace close