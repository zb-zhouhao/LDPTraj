//
// Created by 周昊 on 2023/3/30.
//

#include "ExpMechSolver.h"


ExpMechSolver::ExpMechSolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                             map<int, double> &sens_map, int theta, string parameters_str, map<string, int>& stats_map): Solver(eps, is_out, output_path, trajs_data, sens_map, theta, parameters_str, stats_map) {
    this->threshold = 0.0;
    this->perturbTime = false;
}

void ExpMechSolver::solve() {
    auto it = this->trajsData.begin();
    double consumeTime = 0; double dtw =0;
    this->solveLog.append(Util::getCurTime() + "[INFO] begin randomized trajectories\n");

    int cnt = 0;
    vector<vector<int>> mqeV(Util::queryNum, vector<int>(2, 0));
    while (it != this->trajsData.end()) {
        if (Util::testing && cnt == Util::testCnt) {break;}
        auto start = chrono::high_resolution_clock::now();
        int trajLen = it->second.size();
        double epsilonOverPoint = totalEpsilon / trajLen;
        vector<Triple *> rTraj;
        for (auto &tri : it->second) {
            queue<Triple*> q;
            vector<Triple*> cp;
            vector<double> cost;
            set<long long> vis;
            double sumP = 0;

            q.push(tri);
            vis.insert(tri->X * Util::GRID_SZ + tri->Y);
            long lapT = Util::laplace(tri->time, this->sensMap[it->first], epsilonOverPoint);

            double sens = this->sensMap[it->first];
            int outCnt = 0;
            while (!q.empty()) {
                Triple *tr = q.front();
                q.pop();
                outCnt++;

                int x = tr->X; int y = tr->Y; long time = tr->time;
                double phyDis = Triple::getPhysicDis(tri, x, y);
                double c = exp(epsilonOverPoint * -phyDis / (2 * sens));
                //截断
                if (cost.size() > 0) {
                    double relativeCostRatio = c / cost[0];
                    if (relativeCostRatio < this->threshold) {
                        continue;
                    }
                }
                cost.emplace_back(c);
                sumP += c;
                cp.emplace_back(tr);
                //add candidate
                for (int i = 0; i < 8; i++) {
                    int newX = Util::dirArr[i][0] + x;
                    int newY = Util::dirArr[i][1] + y;
                    if (vis.find(newX * Util::GRID_SZ + newY) != vis.end()) {
                        continue;
                    }

                    if (newX >= 0 && newX < Util::GRID_SZ && newY >= 0 && newY < Util::GRID_SZ) {
                        long newT = time;
                        if (this->perturbTime) {
                            newT = lapT;
                        }
                        auto *tmpTri = new Triple(newT, newX, newY);
                        q.push(tmpTri);
                        vis.insert(newX * Util::GRID_SZ + newY);
                    }

                }
            }
            int retId = Util::chooseFromExpMech(sumP, cp, cost);
            rTraj.emplace_back(cp[retId]);
        }
        if (this->perturbTime) {
            Util::tidyTopCVecTime(rTraj);
        }
        auto end = chrono::high_resolution_clock::now();
        consumeTime += (double)(chrono::duration_cast<chrono::milliseconds>(end - start).count());
//        double dtw_single =  Util::getDTWMetricSingle(it->second, rTraj);
//        cout << "id:" << it->first << " dtw_single: " << dtw_single << endl;
//        dtw += dtw_single;
//        vector<vector<int>> rt;
//        Util::rangeQuery(rTraj, it->second, rt);
//        for (int ii = 0; ii < Util::queryNum; ii++) {
//            mqeV[ii][0] += rt[ii][0];
//            mqeV[ii][1] += rt[ii][1];
//        }
        saveTraj(rTraj, it->first);
        it++;
        cnt++;

    }
    this->solveLog.append(Util::getCurTime() + "[INFO] 算法运行时间: " + to_string(consumeTime / cnt) + " ms\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] finish randomized trajectories\n");
//    this->solveLog.append(Util::getCurTime() + "[INFO] start DTW algorithm to metric similarity between trajectories\n");
//////    cout << cnt << endl;
//    this->solveLog.append(Util::getCurTime() + "[INFO] finish DTW algorithm, DTW value is :" + to_string(dtw / cnt) + "\n");
//    this->solveLog.append(Util::getCurTime() + "[INFO] start MQE algorithm to metric similarity between trajectories\n");
//    this->solveLog.append(Util::getCurTime() + "[INFO] finish MQE algorithm, MQE value is :" + to_string(Util::calMQE(mqeV)) + "\n");
}