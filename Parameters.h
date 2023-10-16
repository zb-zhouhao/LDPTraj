//
// Created by 周昊 on 2023/3/8.
//

#ifndef NNGRAM_PARAMETERS_H
#define NNGRAM_PARAMETERS_H
#include "util.h"
#include "io/FileReader.h"

class Parameters {
public:
    string datasetPath;
    string outputPath;
    int dataTotal;
    int numIndex;
    bool isOut;
    double epsilon;
    vector<double> epsilonVec;
    int roundTotal;
    int theta;
    vector<int> thetaVec;
    double numOfKeyPointRate;
    int runType;

    Parameters();
    string readParametersFromFile(const string &filename);
    string printOut(const map<string, string> &key2value);
};
#endif //NNGRAM_PARAMETERS_H
