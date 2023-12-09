//
// Created by 周昊 on 2023/3/14.
//

#ifndef NNGRAM_FILEWRITER_H
#define NNGRAM_FILEWRITER_H
#include "../util.h"
using namespace std;

class FileWriter {
public:
    static bool outputTrip(const string& filename, vector<Triple *> &traj, map<string, int>& stats_map){
        ofstream outfile;
        outfile.open(filename);
        if (outfile.is_open()) {
            stringstream ss;
            for (auto &tri : traj) {
                ss << to_string(tri->X * stats_map["x_width"] + stats_map["min_x"]) << ", " << to_string(tri->Y * stats_map["y_width"] + stats_map["min_y"]) << ", " << Util::Stamp2Time(tri->time) << "\n";
            }
            if(ss.good()){
                outfile << ss.str();
            }
            ss.clear();
            outfile.close();
            return true;
        }
        return false;
    }

    static bool outputSensMap(const string& filename, map<int, double> &sens_map){
        ofstream outfile;
        outfile.open(filename);
        if (outfile.is_open()) {
            for (const auto& entry : sens_map) {
                outfile << entry.first << ", " << entry.second << std::endl;
            }
            outfile.close(); // 关闭文件
            return true;
        }
        return false;
    }
};


#endif //NNGRAM_FILEWRITER_H
