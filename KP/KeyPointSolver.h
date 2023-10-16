//
// Created by 周昊 on 2023/4/4.
//

#ifndef NNGRAM_KEYPOINTSOLVER_H
#define NNGRAM_KEYPOINTSOLVER_H
#include "../util.h"
#include "../io/FileWriter.h"
#include "../Solver.h"

using namespace std;

class KeyPointSolver : public Solver {
    double numOfKeyPointRate;
    map<int, vector<double>> keyMap;
    map<int, vector<int>> topCKeyPointsMap;
    double threshold;
    bool perturbTime;
    int sampleRate;
    int numOfKeyPoint;
    bool useExpMech;

public:
    KeyPointSolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                   map<int, double> &sens_map, int theta, bool use_exp_mech = true);
    void calKey(vector<Triple *> &traj, vector<double> &keyVec);
    void chooseTopCPoints(vector<int> &topCVec, vector<Triple *> &traj, vector<double> &keyVec);
    void randomizeOneKeyPoint(int keyId, vector<Triple *> &kpRet, vector<Triple *> &traj);
    void randomizeKeyPointsSeq(vector<int> &topCVec, vector<Triple *> &traj, vector<Triple *> &kpRet);
    void genTrajectory(vector<Triple *> &rTraj, vector<Triple *> &kpRet, int sample_rate);
    void solve() override;
    void setNumOfPoint(int traj_len);
    void expSolve(vector<int> &topCVec, vector<Triple *> &traj, vector<Triple *> &kpRet);
    void geoISolve(vector<int> &topCVec, vector<Triple *> &traj, vector<Triple *> &kpRet);
};
#endif //NNGRAM_KEYPOINTSOLVER_H
