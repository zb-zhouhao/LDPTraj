//
// Created by 周昊 on 2023/3/8.
//

#include "NGramSolver.h"

NGramSolver::NGramSolver(double eps, bool is_out, string output_path, map<int, vector<Triple *>> &trajs_data,
                         map<int, double> &sens_map, int theta, string parameters_str, map<string, int>& stats_map) : Solver(eps, is_out, output_path, trajs_data, sens_map, theta, parameters_str, stats_map) {
    initPointBigrams();
}

void NGramSolver::clear() {
    this->timeSeries.clear();
    this->trajs.clear();
    this->trajBigrams.clear();
    this->feasibleBigrams.clear();
}

void NGramSolver::TestInit() {
//    this->timeSeries = {1, 3, 4, 6, 7, 9, 11, 12, 13, 14};
    this->timeSeries = {1, 3, 4, 6, 7};
    this->trajs.emplace_back(0, 1);
    this->trajs.emplace_back(2, 0);
    this->trajs.emplace_back(3, 3);
    this->trajs.emplace_back(4, 3);
    this->trajs.emplace_back(3, 4);

//    this->trajs.emplace_back(5, 7);
//    this->trajs.emplace_back(6, 8);
//    this->trajs.emplace_back(7, 10);
//    this->trajs.emplace_back(9, 13);
//    this->trajs.emplace_back(14, 14);
}

void NGramSolver::initPointBigrams() {
    for (int i = Util::xBegin; i < Util::xEnd; i++) {
        for (int j = Util::yBegin; j < Util::yEnd; j++) {
            this->points.push_back(Util::Point(i, j));
        }
    }
    for (int i = 0; i < this->points.size(); i++) {
        assert(this->points[i].index == i);
    }
    for (auto p1: this->points) {
        for (auto p2: this->points) {
            this->bigrams.push_back(Util::Bigram(p1, p2));
        }
    }
    for (int i = 0; i < this->bigrams.size(); i++) {
        assert(this->bigrams[i].index == i);
    }
};

void NGramSolver::initFeasibleBigrams(int i) {
    double interval = (this->timeSeries[i + 1] - this->timeSeries[i]) / 60.0;
    double thres = (double) interval * this->theta * Util::SCAL_FAC;
    for (auto b: this->bigrams) {
        if (b.getInternalDis() < thres) {
            this->feasibleBigrams[i].push_back(b);
        }
    }
    this->solveLog += "feasible bigrams for " + to_string(i) + " is " + to_string(feasibleBigrams[i].size()) + "\n";
    cout << "feasible bigrams for " << i << " is " << feasibleBigrams[i].size() << endl;
}

int NGramSolver::expMechanismChoose(Util::Bigram ori) {
    default_random_engine randomEngine(time(nullptr));
    double wSum = 0;
    for (auto &b: this->bigrams) {
        double d = ori.getDis(b);
        b.weight = 1.0 / pow(2, d / this->ALPHA);
        wSum += b.weight;
    }
    uniform_real_distribution<double> rand01(0, 1);
    double thres = rand01(randomEngine) * wSum;
    wSum = 0;
    for (int i = 0; i < this->bigrams.size(); i++) {
        wSum += this->bigrams[i].weight;
        if (wSum > thres) {
            return i;
        }
    }
    return this->bigrams.size() - 1;
}

void NGramSolver::initTrajBigrams(int i) {
    Util::Point p1 = this->trajs[i], p2 = this->trajs[i + 1];
    Util::Bigram b(p1, p2);
    int index = expMechanismChoose(b);
    this->solveLog += "choose " + to_string(index) + " for " + to_string(i) + "\n";
    cout << "choose " << index << " for " << i << endl;
    this->trajBigrams[i] = this->bigrams[index];
    this->solveLog += "(" + to_string(this->trajBigrams[i].p1.x) + ", " + to_string(this->trajBigrams[i].p1.y) + "), (" \
 + to_string(this->trajBigrams[i].p2.x) + ", " + to_string(this->trajBigrams[i].p2.y) + ")\n";
    cout << "(" << this->trajBigrams[i].p1.x << ", " << this->trajBigrams[i].p1.y << "), (" <<
         this->trajBigrams[i].p2.x << ", " << this->trajBigrams[i].p2.y << ")" << endl;
}

int NGramSolver::connect(Util::Bigram &b1, Util::Bigram &b2) {
    if (b1.p2.index == b2.p1.index) {
        return 2;
    }
    return 1;
}

int NGramSolver::countConnectConstrainNum() {
    int res = 0;
    for (int i = 1; i < this->feasibleBigrams.size(); i++) {
        res += this->feasibleBigrams[i - 1].size() * this->feasibleBigrams[i].size();
    }
    return res;
}

double NGramSolver::getError(int t, Util::Bigram &b) {
    double res = 0;
    if (t) {
        res += Util::pointDis(this->trajBigrams[t - 1].p2, b.p1);
    }

    res += b.getDis(this->trajBigrams[t]);
    res += Util::pointDis(this->trajBigrams[t + 1].p1, b.p2);
    return res;
}

vector<int> NGramSolver::solveWithPruning() {
    int tNum = this->trajs.size();
    unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
    vector<vector<MPVariable *>> vars(tNum - 1, vector<MPVariable *>(0));
    for (int i = 0; i < this->feasibleBigrams.size(); i++) {
        vector<Util::Bigram> &b = this->feasibleBigrams[i];
        vars[i] = vector<MPVariable *>(b.size());
        for (int w = 0; w < b.size(); w++) {
            vars[i][w] = solver->MakeIntVar(0, 1, "X" + to_string(i) + "@" + to_string(w));
        }
    }
    solveLog += "Number of variables = " + to_string(solver->NumVariables()) + "\n";
    cout << "Number of variables = " << solver->NumVariables() << endl;

    vector<MPConstraint *> ConstraintsForSum(tNum - 1);
    for (int i = 0; i < tNum - 1; i++) {
        ConstraintsForSum[i] = solver->MakeRowConstraint(1, 1, "sum@" + to_string(i + 1));
        for (int j = 0; j < vars[i].size(); j++) {
            ConstraintsForSum[i]->SetCoefficient(vars[i][j], 1);
        }
    }
    this->solveLog += "Number of constraints = " + to_string(solver->NumConstraints()) + "\n";
    cout << "Number of constraints = " << solver->NumConstraints() << endl;

    vector<MPConstraint *> ConstraintsForConnect(countConnectConstrainNum());
    int cur = 0;
    for (int i = 0; i < tNum - 2; i++) {
        for (int w1 = 0; w1 < this->feasibleBigrams[i].size(); w1++) {
            for (int w2 = 0; w2 < this->feasibleBigrams[i + 1].size(); w2++) {
                int v = connect(this->bigrams[w1], this->bigrams[w2]);
                ConstraintsForConnect[cur] = solver->MakeRowConstraint(0, v, "conn@" +
                                                                             to_string(i) + "@" + to_string(w1) + "@" +
                                                                             to_string(w2));
                ConstraintsForConnect[cur]->SetCoefficient(vars[i][w1], 1);
                ConstraintsForConnect[cur]->SetCoefficient(vars[i + 1][w2], 1);
                cur++;
            }
        }
    }
    this->solveLog += "Number of constraints = " + to_string(solver->NumConstraints()) + "\n";
    cout << "Number of constraints = " << solver->NumConstraints() << endl;

    MPObjective *objective = solver->MutableObjective();
    for (int i = 0; i < tNum - 1; i++) {
        for (int w = 0; w < vars[i].size(); w++) {
            double err = getError(i, this->bigrams[w]);
            objective->SetCoefficient(vars[i][w], err);
        }
    }
    objective->SetMinimization();

    this->solveLog += "begin to solve\n";
    cout << "begin to solve\n";
    const MPSolver::ResultStatus result_status = solver->Solve();
//    this->solveLog += to_string(objective->Value()) + "\n";
    this->solveLog += "solve complete\n";
    cout << "solve complete\n";

    vector<int> ret(tNum);
    Util::Bigram bTmp;
    for (int i = 0; i < tNum - 1; i++) {
        for (int w = 0; w < this->feasibleBigrams[i].size(); w++) {
            if (vars[i][w]->solution_value()) {
                bTmp = this->bigrams[feasibleBigrams[i][w].index];
                ret[i] = bTmp.p1.index;
            }
        }
    }
    ret[tNum - 1] = bTmp.p2.index;
    return ret;
}

void NGramSolver::solve() {
    auto it = this->trajsData.begin();
    double consumeTime = 0;

    int cnt = 0;
    while (it != this->trajsData.end()) {
        if (Util::testing && cnt == Util::testCnt) { break; }

        auto start = chrono::high_resolution_clock::now();
        vector<Triple *> rTraj;
        vector<Triple *> oriTraj;
        clear();

        int innerCnt = 0;
        for (auto tri: it->second) {
            if (innerCnt >= 10) {
                break;
            }
            this->timeSeries.emplace_back(tri->time);
            this->trajs.emplace_back(tri->X, tri->Y);
            oriTraj.emplace_back(tri);
            innerCnt++;
        }
//        TestInit();

        this->feasibleBigrams = vector<vector<Util::Bigram>>(this->timeSeries.size() - 1);
        this->trajBigrams = vector<Util::Bigram>(this->timeSeries.size() - 1);
        for (int i = 0; i < this->timeSeries.size() - 1; i++) {
            initFeasibleBigrams(i);
            initTrajBigrams(i);
        }
        auto ans = solveWithPruning();
        int idx = 0;
        for (auto i: ans) {
            this->solveLog += to_string(this->points[i].x) + " " + to_string(this->points[i].y) + "\n";
            cout << this->points[i].x << " " << this->points[i].y << endl;
            auto *tmpTri = new Triple(this->timeSeries[idx++], this->points[i].x, this->points[i].y);
            rTraj.emplace_back(tmpTri);
        }

        auto end = chrono::high_resolution_clock::now();
        consumeTime += (double) (chrono::duration_cast<chrono::milliseconds>(end - start).count());

        saveTraj(rTraj, it->first);
        it++;
        cnt++;
        break;
    }


    this->solveLog.append(Util::getCurTime() + "[INFO] 算法运行时间: " + to_string(consumeTime / cnt) + " ms\n");
    cout << "算法运行时间: " + to_string(consumeTime / cnt) + " ms\n";
    this->solveLog.append(Util::getCurTime() + "[INFO] finish randomized trajectories\n");
}
