//
// Created by 周昊 on 2023/8/3.
//

#ifndef NNGRAM_SOLVER_H
#define NNGRAM_SOLVER_H
#include "util.h"
#include "io/FileWriter.h"

using namespace std;

class Solver {
protected:
    double totalEpsilon;
    int theta;
    string solveLog;
    bool isOut;
    string methodPath;
    map<int, vector<Triple*>> trajsData;
    map<int, double> sensMap;
    string parametersStr;
    map<string, int> statsMap;

public:
    Solver(double eps, bool is_out, string method_path, map<int, vector<Triple *>> &trajs_data,
           map<int, double> &sens_map, int theta, string parameters_str, map<string, int>& stats_map) {
        this->totalEpsilon = eps;
        this->solveLog = "";
        this->isOut = is_out;
        this->methodPath = method_path;
        this->trajsData = trajs_data;
        this->sensMap = sens_map;
        this->theta = theta;
        this->parametersStr = parameters_str;
        this->statsMap = stats_map;
    };

    string getSolveLog() { return this->solveLog; };
    virtual void solve() {};
    void clearLog() { this->solveLog.clear(); };

    void saveTraj(vector<Triple *> &Traj, int userId) {
        if (!this->isOut) {
            return;
        }
        string dirpath = this->methodPath +  "/syn/" + this->parametersStr;
        if (Util::CreateDirectoryIfNotExists(dirpath)) {
            string filename = this->methodPath + "/syn/" + this->parametersStr + "/" + to_string(userId) + ".txt";
            bool status = FileWriter::outputTrip(filename, Traj, this->statsMap);
            if (!status) {
                cout << "[ERROR] save perturbed traj data error" << endl;
                exit(0);
            }
        }
    }
};

#endif //NNGRAM_SOLVER_H
