//
// Created by 周昊 on 2023/4/4.
//

#include "KeyPointSolver.h"

KeyPointSolver::KeyPointSolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                               map<int, double> &sens_map, int theta, bool use_exp_mech): Solver(eps, is_out, output_path, trajs_data,
                                                                               sens_map, theta) {
    this->numOfKeyPointRate = Util::KEY_RATE;
    this->threshold = Util::THRESHOLD;
    this->perturbTime = false;
    this->sampleRate = Util::SAMPLE_RATE;
    this->numOfKeyPoint = 0;
    this->useExpMech = use_exp_mech;
}

void KeyPointSolver::calKey(vector<Triple *> &traj, vector<double> &keyVec) {
    keyVec[0] = 1;
    for (int i = 1; i < traj.size() - 1; i++) {
        keyVec[i] = pow(Util::getKeyOfOnePoint(traj[i - 1], traj[i], traj[i + 1]), 3);
    }
    keyVec[keyVec.size() - 1] = 1;
}

void KeyPointSolver::chooseTopCPoints(vector<int> &topCVec, vector<Triple *> &traj, vector<double> &keyVec) {
    vector<double> cost;
    set<int> cond;
    double sumExp = 0;
    for (int i = 0; i < traj.size(); i++) {
        cond.insert(i);
        double p = exp(((totalEpsilon / 2) * keyVec[i]) / (2 * numOfKeyPoint));
        cost.push_back(p);
        sumExp += p;
    }

    for (int i = 0; i < numOfKeyPoint; i++) {
        //sample from exp distribution
        int idx = Util::chooseFromExpMech(sumExp, cond, cost);
        //update
        cond.erase(idx);
        topCVec.emplace_back(idx);
        sumExp -= cost[idx];
    }
}

void KeyPointSolver::randomizeOneKeyPoint(int keyId, vector<Triple *> &kpRet, vector<Triple *> &traj) {
    vector<Triple *> cp;
    vector<double> cost;
    double sumP = 0;
    long lapT = Util::laplace(traj[keyId]->time, sensMap[keyId], totalEpsilon / (2.0 * numOfKeyPoint));;

    int dis = - (4 * sensMap[keyId] * numOfKeyPoint * log(threshold)) / totalEpsilon;
    int left = Util::maxx(traj[keyId]->X - dis, 0);
    int right = traj[keyId]->X + dis;
    int up = traj[keyId]->Y + dis;
    int bot = Util::maxx(traj[keyId]->Y - dis, 0);
    for (int i = bot; i <= up; i++) {
        for (int j = left; j <= right; j++) {
            long newT = traj[keyId]->time;
            if (perturbTime) {
                newT = lapT;
            }

            auto *tmpTri = new Triple(newT, j, i);
            double r = Triple::getPhysicDis(traj[keyId], tmpTri);
            double score = exp(-(totalEpsilon * r) / (4 * numOfKeyPoint * sensMap[keyId]));
            cp.emplace_back(tmpTri);
            cost.emplace_back(score);
            sumP += score;
        }
    }
    int retId = Util::chooseFromExpMech(sumP, cp, cost);
    kpRet.emplace_back(cp[retId]);
}

void KeyPointSolver::randomizeKeyPointsSeq(vector<int> &topCVec, vector<Triple *> &traj, vector<Triple *> &kpRet) {
    bool fst = true;
    int preId = -1;

    for (auto kid : topCVec) {
        if (fst) {
            randomizeOneKeyPoint(kid, kpRet, traj);
            preId = kid;
            fst = false;
        } else {
            vector<Triple *> cp;
            set<long long> vis;
            vector<double> cost;
            double sumP = 0;
            long lapT = Util::laplace(traj[kid]->time, sensMap[kid], totalEpsilon / (2.0 * numOfKeyPoint));

            int dis = - (4 * sensMap[kid] * numOfKeyPoint * log(threshold)) / totalEpsilon;
            int left = Util::maxx(traj[kid]->X - dis, 0);
            int right = traj[kid]->X + dis;
            int up = traj[kid]->Y + dis;
            int bot = Util::maxx(traj[kid]->Y - dis, 0);
            for (int i = bot; i <= up; i++) {
                for (int j = left; j <= right; j++) {
                    long newT = traj[kid]->time;
                    if (perturbTime) {
                        newT = lapT;
                    }

                    auto *tmpTri = new Triple(newT, j, i);
                    double r = Triple::getPhysicDis(traj[kid], tmpTri);
                    //可达性约束
                    double mins_cost = (traj[kid]->time - traj[preId]->time) / 60.0;
                    if (r > this->theta * Util::SCAL_FAC * mins_cost) {
                        continue;
                    }
                    double score = exp(-(totalEpsilon * r) / (4 * numOfKeyPoint * sensMap[kid]));
                    cp.emplace_back(tmpTri);
                    cost.emplace_back(score);
                    sumP += score;
                }
            }

            preId = kid;
            if (cp.size() == 0) {
                kpRet.emplace_back(traj[kid]);
            } else {
                int retId = Util::chooseFromExpMech(sumP, cp, cost);
                kpRet.emplace_back(cp[retId]);
            }
        }
    }

}

void KeyPointSolver::genTrajectory(vector<Triple *> &rTraj, vector<Triple *> &kpRet, int sample_rate) {
    bool insert = true;
//    bool insert = false;
    for (int k = 0; k < numOfKeyPoint; k++) {
        rTraj.emplace_back(kpRet[k]);
        if (!insert) {
            continue;
        }
        long tim = 1;
        while (k < numOfKeyPoint - 1 && (kpRet[k]->time + tim * sample_rate * 60.0) < kpRet[k + 1]->time) {
            double mins_cost = (kpRet[k + 1]->time - kpRet[k]->time) / 60.0;
            double newX = kpRet[k]->X +
                          (tim * sample_rate) * ((kpRet[k + 1]->X - kpRet[k]->X) / mins_cost);
            double newY = kpRet[k]->Y +
                          (tim * sample_rate) * ((kpRet[k + 1]->Y - kpRet[k]->Y)  / mins_cost);
            long newT = kpRet[k]->time + tim * sample_rate * 60.0;
            auto *tri = new Triple(newT, newX, newY);
            rTraj.emplace_back(tri);
            tim = tim + 1;
        }
    }
}

void KeyPointSolver::solve() {
    auto it = this->trajsData.begin();
    double consumeTime = 0;
    double dtw = 0;
    this->solveLog.append(Util::getCurTime() + "[INFO] begin randomized trajectories\n");

    int cnt = 0;
    vector<vector<int>> mqeV(Util::queryNum, vector<int>(2, 0));
    while (it != this->trajsData.end()) {
        if (Util::testing && cnt == Util::testCnt) { break; }

        auto start = chrono::high_resolution_clock::now();
        //test 轨迹长度为关键点数目
        int trajLen = it->second.size();
        assert(trajLen >= this->numOfKeyPoint);
        this->setNumOfPoint(trajLen);
        //step 1:cal key of each point
        vector<double> kv(trajLen, 0);
        calKey(it->second, kv);
        this->keyMap[it->first] = kv;

        //step 2:choose top c points and sort by time
        //privacy_budget: totalEpsilon / 2
        vector<int> topCVec; //save id
        vector<Triple *> shuffleTopCVec;    //save tri

        chooseTopCPoints(topCVec, it->second, keyMap[it->first]);
        sort(topCVec.begin(), topCVec.end());
        topCKeyPointsMap[it->first] = topCVec;
        assert(topCVec.size() == numOfKeyPoint);

        //step3: randomize top c key point
        //privacy_budget: totalEpsilon / 2
//        randomizeKeyPointsSeq(topCVec, it->second, shuffleTopCVec);
        if (this->useExpMech) {
            expSolve(topCVec, it->second, shuffleTopCVec);
        } else {
            geoISolve(topCVec, it->second, shuffleTopCVec);
        }
        assert(shuffleTopCVec.size() == numOfKeyPoint);

        //sort time
        if (this->perturbTime) {
            Util::tidyTopCVecTime(shuffleTopCVec);
        }

        //step4: traj reconstruct
        vector<Triple *> rTraj; //save gen traj
        genTrajectory(rTraj, shuffleTopCVec, this->sampleRate);

        auto end = chrono::high_resolution_clock::now();
        consumeTime += (double) (chrono::duration_cast<chrono::milliseconds>(end - start).count());
        dtw += Util::getDTWMetricSingle(it->second, rTraj);
        vector<vector<int>> rt;
        Util::rangeQuery(rTraj, it->second, rt);
        for (int ii = 0; ii < Util::queryNum; ii++) {
            mqeV[ii][0] += rt[ii][0];
            mqeV[ii][1] += rt[ii][1];
        }
        saveTraj(rTraj, it->first);
        it++;
        cnt++;
    }
    this->solveLog.append(Util::getCurTime() + "[INFO] 算法运行时间: " + to_string(consumeTime / cnt) + " ms\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] finish randomized trajectories\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] start DTW algorithm to metric similarity between trajectories\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] finish DTW algorithm, DTW value is :" + to_string(dtw / cnt) + "\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] start MQE algorithm to metric similarity between trajectories\n");
    this->solveLog.append(Util::getCurTime() + "[INFO] finish MQE algorithm, MQE value is :" + to_string(Util::calMQE(mqeV)) + "\n");
}

void KeyPointSolver::setNumOfPoint(int traj_len) {
    this->numOfKeyPoint = traj_len * numOfKeyPointRate;
}

void KeyPointSolver::expSolve(vector<int> &topCVec, vector<Triple *> &traj, vector<Triple *> &kpRet) {
    randomizeKeyPointsSeq(topCVec, traj, kpRet);
}

void KeyPointSolver::geoISolve(vector<int> &topCVec, vector<Triple *> &traj, vector<Triple *> &kpRet) {
    double epsilonOverPoint = totalEpsilon / numOfKeyPoint;
    for (auto &tri_id : topCVec) {
        uniform_real_distribution<double> rand02PI(0, 2 * M_PI);
        uniform_real_distribution<double> rand01(0, 1);
        double alpha = rand02PI(randomEngine);
        double p = rand01(randomEngine);
        double r = -1.0 / epsilonOverPoint * (boost::math::lambert_wm1((p - 1) / exp(1)) + 1);
        long newTime = traj[tri_id]->time;
        if (this->perturbTime) {
            newTime = Util::laplace(traj[tri_id]->time, this->sensMap[tri_id], epsilonOverPoint);
        }
        auto *rPoint = new Triple(newTime, traj[tri_id]->X + r * cos(alpha), traj[tri_id]->Y + r * sin(alpha));
        kpRet.emplace_back(rPoint);
    }
}

