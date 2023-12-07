//
// Created by 周昊 on 2023/3/8.
//

#ifndef NNGRAM_FILEREADER_H
#define NNGRAM_FILEREADER_H

#include "../util.h"

using namespace std;

class FileReader {
    ifstream fin;

public:
    explicit FileReader(const string &filename, const string &mode = "open");
    void close();
    static void split(const string &str, vector<string> &tokens, const string &delim);
    void readParameterFile(map<string, string> &key2value, const string &delim);
//    double readTraj(const string &folderName, const int total, map<int, vector<Triple *>>& trajectories, map<int, vector<int>> &recMap, map<int, double> &sensity_map);
    void readTraj(const string &folderName, map<int, vector<Triple *>>& trajectories);
    void statsTraj(map<int, vector<Triple *>>& trajectories, map<string, int>& stats_map);
    void meshTraj(map<int, vector<Triple *>>& trajectories, map<string, int>& stats_map, map<int, double>& sens_map);
};

#endif //NNGRAM_FILEREADER_H
