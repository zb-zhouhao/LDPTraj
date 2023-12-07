//
// Created by 周昊 on 2023/3/8.
//

#ifndef NNGRAM_NGRAMSOLVER_H
#define NNGRAM_NGRAMSOLVER_H
//#include "../util.h"
#include "../Solver.h"

using namespace std;
using namespace operations_research;



class NGramSolver : public Solver {
    vector<long> timeSeries;
    vector<Util::Point> points;
    vector<Util::Bigram> bigrams;
    vector<Util::Point> trajs;
    vector<Util::Bigram> trajBigrams;
    vector<vector<Util::Bigram>> feasibleBigrams;
    const double ALPHA = 0.1;

public:
    NGramSolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                map<int, double> &sens_map, int theta, string parameters_str, map<string, int>& stats_map);
    void clear();
    void TestInit();
    void initPointBigrams();
    void initFeasibleBigrams(int i);
    int expMechanismChoose(Util::Bigram ori);
    void initTrajBigrams(int i);
    int connect(Util::Bigram& b1, Util::Bigram& b2);
    int countConnectConstrainNum();
    double getError(int t, Util::Bigram& b);
    vector<int> solveWithPruning();
    void solve() override;
};

#endif //NNGRAM_NGRAMSOLVER_H
