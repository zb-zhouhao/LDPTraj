//
// Created by 周昊 on 2023/3/14.
//

#ifndef NNGRAM_GEOISOLVER_H
#define NNGRAM_GEOISOLVER_H
#include "../util.h"
#include "../io/FileWriter.h"
#include "../Solver.h"

using namespace std;

class GeoISolver: public Solver {
    bool perturbTime;
public:
    GeoISolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
               map<int, double> &sens_map, int theta, string parameters_str, map<string, int>& stats_map);
    void solve() override;
};

#endif //NNGRAM_GEOISOLVER_H
