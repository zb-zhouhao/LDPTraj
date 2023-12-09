#include "util.h"
#include "io/FileReader.h"
#include "Parameters.h"
#include "io/LogFile.h"
#include "SolverFactory.h"

string printOut(const map<string, int> &key2value)  {
    string msg = "[dataset stats]\n";
    map<string, string>::iterator itr;
    for (auto &pair: key2value) {
        msg.append(pair.first);
        msg.append(" = ");
        msg.append(to_string(pair.second));
        msg.append("\n");
    }
    msg.append("---------------------------------------------------------------------------\n");
    return msg;
}


int main() {
    //parameters
    auto *param = new Parameters();
    string msg = param->readParametersFromFile("config");
    const string datasetPath = param->datasetPath;
    const string outputPath = param->outputPath;
    const double epsilon = param->epsilon;
    const vector<double> epsilonV = param->epsilonVec;
    const int roundTotal = param->roundTotal;
    const int numIndex = param->numIndex;
    const bool isOut = param->isOut;
    const int theta = param->theta;
    const string epsStr = param->epsStr;
    const vector<int> thetaV = param->thetaVec;
    const vector<int> splitIndex = param->splitIndex;

    // "GEO-I"
//    const string LDPMethod = "NG";
//    const string LDPMethod = "Geo-I";
    const string LDPMethod = "Local";
//    const string LDPMethod = "KP-E";
//    const string LDPMethod = "KP-G";
    delete param; param = nullptr;

    // logfile
    for (int idx : splitIndex) {
        const string parametersStr = "n_" + to_string(idx) + "_eps_" + epsStr + "_the_" + to_string(theta);
//    cout << parametersStr << endl;
        const string logFilePath = LogFile::getLogFileName(outputPath, parametersStr, LDPMethod);
        auto *logger = new LogFile(logFilePath);
        logger->addContent(msg);

        //read trajs data
        map<int, vector<Triple *> > trajectories;
        map<int, double> sensityMap;
        map<string, int> statsMap;
        string traj_path = datasetPath + "/geolife_bj_" + to_string(idx);

        msg.append(Util::getCurTime() + "[INFO] reading trajectory data from: " + traj_path + "\n");
        logger->addContent(msg);


        //load sensMap

        FileReader *freader = new FileReader(traj_path, "exist");
//    freader->readTraj(datasetPath, dataTotal, trajectories, recMap, sensityMap);
        freader->readTraj(traj_path, trajectories);
        freader->statsTraj(trajectories, statsMap);
        freader->meshTraj(trajectories, statsMap, sensityMap);

        msg.append(printOut(statsMap) + "\n");
        logger->addContent(msg);

        //工厂模式
        auto *slvr = SolverFactory::createSolver(LDPMethod, epsilon, isOut, outputPath, trajectories, sensityMap, msg, theta, parametersStr, statsMap);
        for (int rd = 1; rd <= roundTotal; rd++) {
            msg.append( "**********************************************************\n");
            msg.append("[INFO] ROUND: " + to_string(rd) + "\n");
            slvr->solve();
            msg.append(move(slvr->getSolveLog()));
            logger->addContent(msg);
            slvr->clearLog();
        }

        delete slvr; slvr = nullptr;
    }

    return 0;
}

