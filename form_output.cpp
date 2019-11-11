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
                // std::cout << "from " << elem[0] << " to " << elem[1] << " ratio" <<
                // elem[2] << std::endl;
                out.push_back(elem);
                elem.clear();
            }
        }
        return (out);
    }

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
//    std::map<std::string, std::vector<std::vector

} // namespace close