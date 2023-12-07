//
// Created by 周昊 on 2023/3/8.
//

#include "Parameters.h"
Parameters::Parameters() {
    this->datasetPath = "";
    this->outputPath = "";
    this->dataTotal = 0;
    this->numIndex = 0;
    this->epsilon = 0;
    this->isOut = false;
    this->roundTotal = 0;
    this->theta = 0;
    this->runType = 0;
    this->numOfKeyPointRate = 0.0;
    this->epsStr = "";
}

char *getCurrDirectory() {
    char *currDirectory;
    if ((currDirectory = getcwd(nullptr, 0)) == nullptr) {
        perror("getcwd error");
    }
    else {
        printf("Current Working Directory = %s\n", currDirectory);
    }
    return currDirectory;
}

string Parameters::readParametersFromFile(const string &filename) {
    char *currDirectory = getCurrDirectory();

    string fullpath(currDirectory);
    fullpath.append("/");
    fullpath.append(filename);    // fixed

    auto *freader = new FileReader(fullpath);
    map<string, string> key2value;
    freader->readParameterFile(key2value, " = ");

    for (auto &pair: key2value) {
        if (pair.first == "DATASET_PATH") {
            datasetPath = pair.second;
        }
        else if (pair.first == "OUTPUT_PATH") {
            outputPath = pair.second;
        }
        else if (pair.first == "DATA_TOTAL") {
            dataTotal = stoi(pair.second);
        }
        else if (pair.first == "NUM_INDEX") {
            numIndex = stoi(pair.second);
        }
        else if (pair.first == "IS_OUT") {
            isOut = (pair.second == "true");
        }
        else if (pair.first == "EPSILON") {
            epsilon = stod(pair.second);
            epsStr = pair.second;
        }
        else if (pair.first == "EPSILON_V") {
            vector<string> tmp;
            FileReader::split(pair.second, tmp, ",");
            for (const string &v: tmp) {
                epsilonVec.emplace_back(stod(v));
            }
        }
        else if (pair.first == "ROUND_TOTAL") {
            roundTotal = stoi(pair.second);
        }
        else if (pair.first == "THETA") {
            theta = stoi(pair.second);
        }
        else if (pair.first == "THETA_V") {
            vector<string> tmp;
            FileReader::split(pair.second, tmp, ",");
            for (const string &v: tmp) {
                thetaVec.emplace_back(stoi(v));
            }
        }
        else if (pair.first == "NUM_OF_KEY_POINT_RATE") {
            this->numOfKeyPointRate = stod(pair.second);
        }
        else if (pair.first == "RUN_TYPE") {
            runType  = stoi(pair.second);
        }
        else if (pair.first == "SPLIT_INDEX") {
            vector<string> tmp;
            FileReader::split(pair.second, tmp, ",");
            for (const string &v: tmp) {
                splitIndex.emplace_back(stoi(v));
            }
        }
    }
    return printOut(key2value);
}

string Parameters::printOut(const map<string, string> &key2value)  {
    string msg = "[Parameters]\n";
    map<string, string>::iterator itr;
    for (auto &pair: key2value) {
        msg.append(pair.first);
        msg.append(" = ");
        msg.append( pair.second);
        msg.append("\n");
    }
    msg.append("---------------------------------------------------------------------------\n");
    return msg;
}
