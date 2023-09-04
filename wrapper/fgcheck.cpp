//
// Created by Yifan Chen on 2023/9/1.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <vector>

using namespace std;

int main(int argc, char * argv[]) {
    ifstream fg_ifs;
    fg_ifs.open(argv[1]);
    size_t num_factors;
    fg_ifs >> num_factors;

    set<size_t> all_vars;
    multimap<size_t, size_t> prec;
    set<size_t> bodies;

    for (size_t f_id = 0; f_id < num_factors; ++f_id) {
        size_t num_f_vars;
        fg_ifs >> num_f_vars;

        size_t head;
        fg_ifs >> head;
        all_vars.insert(head);

        for (size_t v_id = 1; v_id < num_f_vars; ++v_id) {
            size_t body;
            fg_ifs >> body;
            all_vars.insert(body);
            prec.emplace(head, body);
            bodies.insert(body);
        }
        for (size_t v_id = 0; v_id < num_f_vars; ++v_id) {
            size_t dim;
            fg_ifs >> dim;
        }
        size_t num_col;
        fg_ifs >> num_col;
        for (size_t col_id = 0; col_id < num_col; ++col_id) {
            size_t col;
            double val;
            fg_ifs >> col >> val;
        }
    }
    fg_ifs.close();
    clog << "All variables #: " << all_vars.size() << endl;
    clog << "Terminal variables #: " << all_vars.size() - bodies.size() << endl;
    for (size_t v : all_vars) {
        if (bodies.find(v) == bodies.end()) {
            clog << v << "\t";
        }
    }
    clog << endl;

    ifstream tab_ifs;
    tab_ifs.open(argv[2]);
    string line;
    getline(tab_ifs, line);
    istringstream var_ss(line);
    size_t obs_id;
    vector<size_t> obs;
    while (var_ss >> obs_id) {
        obs.push_back(obs_id);
    }
    clog << "Total observable variables #: " << obs.size() << endl;
    getline(tab_ifs, line);

    set<size_t> conditioned_once;

    while (getline(tab_ifs, line)) {
        vector<size_t> cur_obs;
        istringstream ob_ss(line);
        for (auto ob : obs) {
            string ob_str;
            getline(ob_ss, ob_str, '\t');
            if (!ob_str.empty()) {
                cur_obs.push_back(ob);
            }
        }
        clog << "\tObserved variables #: " << cur_obs.size() << std::endl;
        deque<size_t> worklist{cur_obs.begin(), cur_obs.end()};
        set<size_t> conditioned(cur_obs.begin(), cur_obs.end());
        while (!worklist.empty()) {
            obs_id = worklist.front();
            worklist.pop_front();
            auto prec_range = prec.equal_range(obs_id);
            for (auto it = prec_range.first; it != prec_range.second; ++it) {
                auto body = it->second;
                if (conditioned.find(body) == conditioned.end()) {
                    conditioned.insert(body);
                    worklist.push_back(body);
                }
            }
        }
//        for (auto cond_var: conditioned) {
//            all_vars.erase(cond_var);
//        }

        clog << "\tUnconditioned variables #: " << all_vars.size() - conditioned.size() << std::endl;
        conditioned_once.insert(conditioned.begin(), conditioned.end());
    }

    clog << "Never conditioned variables #: " << all_vars.size() - conditioned_once.size() << std::endl;
    tab_ifs.close();
    return 0;
}