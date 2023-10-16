#include "util.h"
#include "io/FileReader.h"
#include "Parameters.h"
#include "io/LogFile.h"
#include "SolverFactory.h"

int main() {
    //parameters
    auto *param = new Parameters();
    string msg = param->readParametersFromFile("config");
    const string datasetPath = param->datasetPath;
    const string outputPath = param->outputPath;
    const double epsilon = param->epsilon;
    const vector<double> epsilonV = param->epsilonVec;
    const int roundTotal = param->roundTotal;
    const int dataTotal = param->dataTotal;
    const int numIndex = param->numIndex;
    const bool isOut = param->isOut;
    const int theta = param->theta;
    const vector<int> thetaV = param->thetaVec;
    const int runType = param->runType;

    // "GEO-I"
//    const string LDPMethod = "NG";
//    const string LDPMethod = "Geo-I";
//    const string LDPMethod = "Local";
//    const string LDPMethod = "KP-E";
    const string LDPMethod = "KP-G";
    delete param; param = nullptr;

    // logfile
    const string parametersStr = "_RT" + to_string(runType) + "_D" + to_string(dataTotal) + "_IDX" + to_string(numIndex) ;
    const string logFilePath = LogFile::getLogFileName(outputPath, parametersStr, LDPMethod);
    auto *logger = new LogFile(logFilePath);
    logger->addContent(msg);

    //read trajs data
    map<int, vector<Triple *> > trajectories;
    map<int, double> sensityMap;

    //load sensMap

    msg.append(Util::getCurTime() + "[INFO] reading trajectory  data from: " + datasetPath + "\n");
    logger->addContent(msg);

    auto *freader = new FileReader(datasetPath, "exist");
    map<int, vector<int>> recMap;
    double avgL = freader->readTraj(datasetPath, dataTotal, trajectories, recMap, sensityMap);
    for (auto it = recMap.begin(); it != recMap.end(); it++) {
        msg.append(Util::getCurTime() + "[INFO] num of trajectories len in " + to_string(it->first) + "+-(25):" + to_string(it->second.size()) + "\n");
    }
    msg.append(Util::getCurTime() + "[INFO] average len: " + to_string(avgL) + "\n");
    logger->addContent(msg);

    //get traj data from index
    map<int, vector<Triple *> > testTrajectories;
    for (auto idx : recMap[numIndex]) {
        testTrajectories[idx] = trajectories[idx];
    }
    trajectories.clear();
    msg.append(Util::getCurTime() + "[INFO] numIndex: " + to_string(numIndex) + "\n");
    msg.append(Util::getCurTime() + "[INFO] trajs size: " + to_string(testTrajectories.size()) + "\n");
    logger->addContent(msg);

    //save traj data if needed
    for (auto it = testTrajectories.begin(); it != testTrajectories.end(); it++) {
        bool status = FileWriter::outputTrip(outputPath + "ori/" + to_string(it->first) + ".txt", it->second);
        if (!status) {
            msg.append(Util::getCurTime() + "[ERROR] save ori traj data error\n");
            logger->addContent(msg);
            exit(0);
        }
    }

    //工厂模式
    auto *slvr = SolverFactory::createSolver(LDPMethod, epsilon, isOut, outputPath, testTrajectories, sensityMap, msg, theta);
    for (int rd = 1; rd <= roundTotal; rd++) {
        msg.append( "**********************************************************\n");
        msg.append("[INFO] ROUND: " + to_string(rd) + "\n");
        slvr->solve();
        msg.append(move(slvr->getSolveLog()));
        logger->addContent(msg);
        slvr->clearLog();
    }
    delete slvr; slvr = nullptr;
    return 0;
}
