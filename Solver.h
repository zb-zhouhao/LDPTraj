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
    string outputPath;
    map<int, vector<Triple*>> trajsData;
    map<int, double> sensMap;

public:
    Solver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
           map<int, double> &sens_map, int theta) {
        this->totalEpsilon = eps;
        this->solveLog = "";
        this->isOut = is_out;
        this->outputPath = output_path;
        this->trajsData = trajs_data;
        this->sensMap = sens_map;
        this->theta = theta;
    };

    string getSolveLog() { return this->solveLog; };
    virtual void solve() {};
    void clearLog() { this->solveLog.clear(); };
    void saveTraj(vector<Triple *> &Traj, int userId) {
        if (!this->isOut) {
            return;
        }
        string filename = this->outputPath + "per/" + to_string(userId) + ".txt";
        bool status = FileWriter::outputTrip(filename, Traj);
        if (!status) {
            cout << "[ERROR] save perturbed traj data error" << endl;
            exit(0);
        }
    }
};

#endif //NNGRAM_SOLVER_H
