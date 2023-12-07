//
// Created by 周昊 on 2023/8/3.
//

#ifndef NNGRAM_SOLVERFACTORY_H
#define NNGRAM_SOLVERFACTORY_H

#include "Solver.h"
#include "NG/NGramSolver.h"
#include "Local/ExpMechSolver.h"
#include "GI/GeoISolver.h"
#include "KP/KeyPointSolver.h"

class SolverFactory {
public:
    static Solver* createSolver(string method, double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                 map<int, double> &sens_map, string& msg, int theta, string parameters_str, map<string, int>& stats_map

    ) {
        Solver *slvr = nullptr;
        msg.append("-----------------------------------------------------------------------------\n");
        string method_path = output_path.append(method);
        if (method == "NG") {
            msg.append("-----------------------------[ NGram LDP ]----------------------------------\n");

            slvr = new NGramSolver(eps, is_out, method_path, trajs_data, sens_map, theta, parameters_str, stats_map);
        } else if (method == "Geo-I") {
            msg.append("-----------------------------[ Geo-I LDP ]----------------------------------\n");
            slvr = new GeoISolver(eps, is_out, method_path, trajs_data, sens_map, theta, parameters_str, stats_map);
        } else if (method == "Local") {
            msg.append("-----------------------------[ Local LDP ]----------------------------------\n");
            slvr = new ExpMechSolver(eps, is_out, method_path, trajs_data, sens_map, theta, parameters_str, stats_map);
        } else if (method == "KP-E") {
            bool use_exp_mech = true;
            msg.append("-----------------------------[ KP-E  LDP ]----------------------------------\n");
            slvr = new KeyPointSolver(eps, is_out, method_path, trajs_data, sens_map, theta, stats_map, use_exp_mech, parameters_str);
        } else if (method == "KP-G") {
            bool use_exp_mech = false;
            msg.append("-----------------------------[ KP-E  LDP ]----------------------------------\n");
            slvr = new KeyPointSolver(eps, is_out, method_path, trajs_data, sens_map, theta, stats_map, use_exp_mech, parameters_str);
        }
        msg.append("-----------------------------------------------------------------------------\n");
        return slvr;
    }

};

#endif //NNGRAM_SOLVERFACTORY_H
