//
// Created by 周昊 on 2023/3/8.
//

#include "FileReader.h"

// constructor
FileReader::FileReader(const string &filename, const string &mode) {
    if (mode == "open") {
        fin.open(filename, ifstream::in);
        if (!fin.good()) {
            string msg("file not found: ");
            msg.append(Util::getCurTime() + filename);
            throw runtime_error(msg);
        }
    }
    else if(mode == "exist"){   // whether the folder exists
        if (access(filename.c_str(), 0) == -1){
//        if (_access(filename.c_str(), 0) == -1){
            string msg("folder does not exist: ");
            msg.append(Util::getCurTime() + filename);
            throw runtime_error(msg);
        }
    }
}

void FileReader::close() {
    if (fin.is_open()) {
        fin.close();
    }
}

void FileReader::split(const string &str, vector<string> &tokens, const string &delim) {
    tokens.clear();
    auto start = str.find_first_not_of(delim, 0);
    auto position = str.find_first_of(delim, start);
    while (position != string::npos || start != string::npos) {
        tokens.emplace_back(move(str.substr(start, position - start)));        // [start, position)
// next token
        start = str.find_first_not_of(delim, position);
        position = str.find_first_of(delim, start);
    }
}

void FileReader::readParameterFile(map<string, string> &key2value, const string &delim) {
    string line, paramName, paramValue;
    while (fin.good()) {
        getline(fin, line);
        if (!line.empty() && line[0] != '#') {
            vector<string> vec;
            split(line, vec, delim);
            paramName = vec[0];
            paramValue = vec[1];
            key2value[paramName] = paramValue;
        }
    }
    fin.close();
}

long timestr_timestamp(const string &datetime = "1970-01-01 00:00:00") {
    if (datetime.length() < 19) {
        cout << "invalid string - cant convert to timestamp\n";
    }

    struct tm tm{};
    bool discretization = false;
    tm.tm_year = atoi(datetime.substr(0, 4).c_str()) - 1900;
    tm.tm_mon = atoi(datetime.substr(5, 2).c_str()) - 1;
    tm.tm_mday = atoi(datetime.substr(8, 2).c_str());
    tm.tm_hour = atoi(datetime.substr(11, 2).c_str());
    tm.tm_min = atoi(datetime.substr(14, 2).c_str());
    //分钟对齐
    if (!discretization) {
        tm.tm_sec = atoi(datetime.substr(17, 2).c_str());
    } else {
        string sec = "00";
        tm.tm_sec = atoi(sec.c_str());
    }
    return mktime(&tm);
}

double FileReader::readTraj(const string &folderName, const int total, map<int, vector<Triple *>> &trajectories, map<int, vector<int>> &recMap, map<int, double> &sensity_map) {
    // only for geolife_bj
    const int startId = 0;
    const int laIndex = 0, loIndex = 1, tIndex = 2;
//    map<int, int> recMap;
//    int arr[6] = {7, 50, 100, 500, 1000, 5000};
    int arr[5] = {50, 100, 500, 1000, 5000};
//    int arr[1] = {7};/

    int userId = startId;
    while (trajectories.size() < total) {
//        cout << userId << endl;
//        Util::printProgressBar(userId, total);
        long maxX = -1, minX = LONG_MAX, maxY = -1, minY = LONG_MAX;
        string filename(folderName);
        filename.append(to_string(userId) + ".txt");
        ifstream fs(filename);
        if (!fs.good()) {
            userId++;
            continue;
        }

        vector<Triple *> tra;
        string line;
        long preTimestamp = -1;
        while(getline(fs, line)) {

            vector<string> tokens;
            split(line, tokens, ",");

            if (tokens.size() == 2) {
                continue;
            }
            //去除连续点
            long curTimestamp = timestr_timestamp(tokens.at(tIndex));
            if(preTimestamp < 0 || curTimestamp != preTimestamp) {   // noise: some consecutive points have the same timestamp
                preTimestamp = curTimestamp;
//                int x = (int)(stod(tokens[laIndex]) / Util::X_WIDTH);
//                int y = (int)(stod(tokens[loIndex]) / Util::Y_WIDTH);
                int x = (int)(stod(tokens[laIndex]) / Util::UNIT);
                int y = (int)(stod(tokens[loIndex]) / Util::UNIT);
                minX = minX < x ? minX : x;
                maxX = maxX > x ? maxX : x;
                minY = minY < y ? minY : y;
                maxY = maxY > y ? maxY : y;
                auto *triple = new Triple(curTimestamp, x, y);
                tra.emplace_back(triple);
            }
        }
        int sz = tra.size();
        for (auto p : arr) {
            if (sz <= p + 25 && sz >= p - 25) {
                if (!recMap[p].empty()) {
                    recMap[p].push_back(userId);
                } else {
                    vector<int> indexVec = {userId};
                    recMap[p] = indexVec;
                }
            }
        }
        trajectories[userId] = tra;
        //使用rtree求解敏感度，及轨迹两点之间的最大值
        sensity_map[userId] = sqrt((double)(pow(maxX - minX, 2) + pow(maxY - minY, 2)));
        userId++;
    }

    double avgLength = 0;
    auto it = trajectories.begin();
    while (it != trajectories.end()) {
        avgLength += (double)it->second.size() / (double)trajectories.size();
        it++;
    }
//    cout << trajectories.size() << endl;
    return avgLength;
}
