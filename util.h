//
// Created by 周昊 on 2023/3/8.
//

#ifndef NNGRAM_UTIL_H
#define NNGRAM_UTIL_H
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <random>
#include <ctime>
#include <iomanip>
#include <string>
#include <fstream>
#include <random>
#include <map>
#include <queue>
#include <set>
#include "assert.h"
#include "cmath"
#include <sys/stat.h>
#include <filesystem>

#include "ortools/linear_solver/linear_solver.h"
#include "boost/math/special_functions/lambert_w.hpp"
#include "boost/geometry.hpp"
#include "boost/geometry/index/rtree.hpp"

#include "Triple.h"


using namespace std;
namespace fs = std::filesystem;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef std::pair<point, unsigned> value;

static default_random_engine randomEngine(time(nullptr));

inline int gcd(int a,int b) {
    return b>0 ? gcd(b,a%b):a;
}

class Util {
public:
    static constexpr int dirArr[8][2] = {
            {-1, -1},
            {-1, 0},
            {-1, 1},
            {0, -1},
            {0, 1},
            {1, -1},
            {1, 0},
            {1, 1}
    };

    static constexpr bool testing = false;
    static constexpr int testCnt = 10;

    //自适应读入
    static constexpr int GRID_SZ = 30;
    static constexpr int X_WIDTH = 1;
    static constexpr int Y_WIDTH = 1;

    static constexpr int V = 8; //km/h
    static constexpr double VELOCITY = (V / 3.6);


    static constexpr int xBegin = 0, yBegin = 0, xEnd = GRID_SZ, yEnd = GRID_SZ;
    static constexpr int xLength = xEnd - xBegin, yLength = yEnd - yBegin;
    static constexpr int region_size = xLength * yLength;
    /* ldp算法参数设置 */
    static constexpr int UNIT = 1000; //网格单位：1km
    static constexpr double SCAL_FAC = (1.0 / 60);  //km/min
    static constexpr double KEY_RATE = 0.3;
    static constexpr double THRESHOLD = 0.99;
//    static constexpr int SAMPLE_SEED = ;
//    static constexpr int SAMPLE_RATE = 5 * KEY_RATE * ;
//    static constexpr int SAMPLE_RATE = 5;
    static constexpr double SAMPLE_RATE = 0.2; // 2/3 * 0.3 40s
    struct Point {
        int x, y;
        int index;
        Point(double xx, double yy) {
            x = (int)(xx / X_WIDTH);
            y = (int)(yy / Y_WIDTH);
            index = x * xLength + y;
        }
        Point(int x, int y): x(x), y(y) {
            index = x * xLength + y;
        }
        bool operator==(const Point& p) const {
            return index == p.index;
        }
    };

    static double pointDis(const Point& p1, const Point& p2) {
        return sqrt((double)((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)));
    }

    struct Bigram {
        Point p1, p2;
        int index;
        double weight;
        Bigram(): p1({0, 0}), p2({0, 0}), index(0), weight(0.0) {}
        Bigram(Point p1, Point p2): p1(p1), p2(p2), weight(0.0) {
            index = p1.index * region_size + p2.index;
        }
        double getInternalDis() const {
            return pointDis(p1, p2);
        }
        double getDis(const Bigram& b) const {
            return pointDis(p1, b.p1) + pointDis(p2, b.p2);
        }
    };

    static string Stamp2Time(long timestamp) {
        time_t tick = (time_t) timestamp;
        struct tm tm;
        char s[40];
        tm = *localtime(&tick);
        strftime(s, sizeof(s), "%Y-%m-%d %H:%M:%S", &tm);
        string str(s);
        return str;
    }

    static long timestr_timestamp(const char* c_str) {
        string datetime(c_str);
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

    static string getCurTime() {
        std::time_t currentTime = std::time(nullptr);
        std::tm localTime = *std::localtime(&currentTime);

        // Format the time as a string
        std::stringstream ss;
        ss << localTime.tm_year + 1900 << "-"
           << std::setfill('0') << std::setw(2) << localTime.tm_mon + 1 << "-"
           << std::setfill('0') << std::setw(2) << localTime.tm_mday << " "
           << std::setfill('0') << std::setw(2) << localTime.tm_hour << ":"
           << std::setfill('0') << std::setw(2) << localTime.tm_min << ":"
           << std::setfill('0') << std::setw(2) << localTime.tm_sec << " ";

        return ss.str();
    }

    //input is vec
    static double getDotProduct(const Triple* v1, const Triple* v2) {
        return v1->X * v2->X + v1->Y * v2->Y;
    }

    //input is vec
    static double getVecMod(const Triple* v) {
        return sqrt(pow(v->X, 2) + pow(v->Y, 2));
    }



    //tidy
    static void tidy(Triple *tr) {
        int x = tr->X;
        int y = tr->Y;
        int gcdd = gcd(x, y);
        if (gcdd == 0 || gcdd == 1) {
            return;
        }
        tr->X = x / gcdd;
        tr->Y = y / gcdd;
    }


    //get key of a trajectory
    static double getKeyOfOnePoint(Triple *pre, Triple *cur, Triple *post) {
        Triple *v1 = new Triple(cur->time - pre->time, cur->X - pre->X, cur->Y - pre->Y);
        Triple *v2 = new Triple(post->time - cur->time, post->X - cur->X, post->Y - cur->Y);
        tidy(v1); tidy(v2);
        int v1x = v1->X; int v1y = v1->Y;
        int v2x = v2->X; int v2y = v2->Y;
        double dotProduct = getDotProduct(v1, v2);
        double mod1 = getVecMod(v1);
        double mod2 = getVecMod(v2);
//        cout << v1->X << " " << v2->Y << endl;
        delete v1; delete v2;
        //同向
        if (v1x == v2x && v1y == v2y) {
            return 0;
        } else if (mod1 == 0 && mod2 == 0) {
            return 0;
        } else if ((mod1 == 0 && mod2 != 0) || (mod1 != 0 && mod2 == 0)) {
            return 1;
        } else {
            if (pow(dotProduct / (mod1 * mod2), 2) > 1) {
                cout << "ERROR" << endl;
                exit(0);
            }
            return sqrt(1 - pow(dotProduct / (mod1 * mod2), 2));
        }
    }

    static int chooseFromExpMech(double sumExp, set<int> &cond, vector<double> &cost) {
        uniform_real_distribution<double> rand01(0, 1);
        double threshold = rand01(randomEngine) * sumExp;

        double curSum = 0;
        for (auto idx : cond) {
            curSum += cost[idx];
            if (curSum > threshold) {
                return idx;
            }
        }
        return *cond.rbegin();
    }

    static int chooseFromExpMech(double sumExp, vector<Triple*> &cond, vector<double> &cost) {
        uniform_real_distribution<double> rand01(0, 1);
        double threshold = rand01(randomEngine) * sumExp;

        double curSum = 0;
        for (int i = 0; i < cond.size(); i++) {
            curSum += cost[i];
            if (curSum > threshold) {
                return i;
            }
        }

        return cond.size() - 1;
    }

    static long maxx(long l1, long l2) {
        return l1 < l2 ? l2 : l1;
    }

    static long minn(long l1, long l2) {
        return l1 < l2 ? l1 : l2;
    }

    static long laplace(long oriTime, long sensity, double epsilon) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
//        cout << epsilon << endl;
        long scale = sensity / epsilon;
        double u = dis(gen) - 0.5;  // 生成一个介于 -0.5 和 0.5 之间的随机数
        double noise = -scale * std::copysign(1.0, u) * std::log(1 - 2 * std::abs(u));
//        cout << noise / 60.0 << endl;
        return (long)(oriTime + (noise / 60.0));
    }

    static void tidyTopCVecTime(vector<Triple *> &traj) {
        vector<long> times; int n = traj.size();
        for (auto *tr : traj) {
            times.push_back(tr->time);
        }
        sort(times.begin(), times.end());
        for (int i = 0; i < n; i++) {
            traj[i]->time = times[i];
        }
    }

    static void printProgressBar(int progress, int total, int barWidth = 70) {
        float percentage = static_cast<float>(progress) / total;
        int progressWidth = static_cast<int>(percentage * barWidth);
        cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < progressWidth)
                cout << "=";
            else
                cout << " ";
        }
        cout << "] " << static_cast<int>(percentage * 100.0) << "%\r";
        cout.flush();
    }

    static bool CreateDirectoryIfNotExists(const std::string& directoryPath) {
        int status = mkdir(directoryPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        if (status == 0) {
            return true; // 目录创建成功
        } else if (errno == EEXIST) {
            return true; // 目录已存在
        } else {
            return false; // 目录创建失败
        }
    }
};

#endif //NNGRAM_UTIL_H
