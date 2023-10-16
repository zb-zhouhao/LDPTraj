//
// Created by 周昊 on 2023/3/14.
//

#include "GeoISolver.h"


GeoISolver::GeoISolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                       map<int, double> &sens_map, int theta) : Solver(eps, is_out, output_path, trajs_data, sens_map, theta
) {
    this->perturbTime = false;
}

void GeoISolver::solve() {
    auto it = this->trajsData.begin();
    double consumeTime = 0;
    double dtw = 0;

    this->solveLog.append(Util::getCurTime() + "[INFO] begin randomized trajectories\n");

    int cnt = 0;
    while (it != this->trajsData.end()) {
        if (Util::testing && cnt == Util::testCnt) { break; }
        auto start = chrono::high_resolution_clock::now();
        int trajLen = it->second.size();
        double epsilonOverPoint = totalEpsilon / trajLen;
        vector<Triple *> rTraj;
        for (auto &tri : it->second) {
            //add plane laplace noise
//            default_random_engine randomEngine(time(nullptr));
            //geoi
            uniform_real_distribution<double> rand02PI(0, 2 * M_PI);
            uniform_real_distribution<double> rand01(0, 1);
            double theta = rand02PI(randomEngine);
            double p = rand01(randomEngine);
            double r = -1.0 / epsilonOverPoint * (boost::math::lambert_wm1((p - 1) / exp(1)) + 1);
            long newTime = tri->time;
            if (this->perturbTime) {
                newTime = Util::laplace(tri->time, this->sensMap[it->first], epsilonOverPoint);
            }
            auto *rPoint = new Triple(newTime, tri->X + r * cos(theta), tri->Y + r * sin(theta));
            rTraj.emplace_back(rPoint);
        }

        auto end = chrono::high_resolution_clock::now();
        consumeTime += (double) (chrono::duration_cast<chrono::milliseconds>(end - start).count());
        dtw += Util::getDTWMetricSingle(it->second, rTraj);
        cout << "gen user " << it->first << " randomized traj" << endl;
        saveTraj(rTraj, it->first);
        it++;
        cnt++;

    }

    this->solveLog.append(Util::getCurTime() + "[INFO] 算法运行时间: " + to_string(consumeTime / cnt) + " ms\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] finish randomized trajectories\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] start DTW algorithm to metric similarity between trajectories\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] finish DTW algorithm, DTW value is :" + to_string(dtw / cnt) + "\n");
}


