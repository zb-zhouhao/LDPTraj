//
// Created by 周昊 on 2023/3/30.
//

#ifndef NNGRAM_EXPMECHSOLVER_H
#define NNGRAM_EXPMECHSOLVER_H
#include "../util.h"
#include "../io/FileWriter.h"
#include "../Solver.h"

using namespace std;
using namespace operations_research;

class ExpMechSolver : public Solver {
    double threshold;
    bool perturbTime;

public:
    ExpMechSolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                  map<int, double> &sens_map, int theta);
    void solve() override;
};
#endif //NNGRAM_EXPMECHSOLVER_H
