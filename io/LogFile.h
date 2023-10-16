//
// Created by 周昊 on 2023/3/8.
//

#ifndef NNGRAM_LOGFILE_H
#define NNGRAM_LOGFILE_H

#include "../util.h"

using namespace std;

class LogFile {
    string logFilePath;

public:
    static bool testOrCreateDir(const string &pathName);

    explicit LogFile(string filePath);

    void addContent(string &content, bool screenPrint = true);

    static string getLogFileName(const string &outputFolder, const string &parametersStr, \
                                    const string &LDPMethod);

};

#endif //NNGRAM_LOGFILE_H
