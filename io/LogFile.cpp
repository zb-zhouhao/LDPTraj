//
// Created by 周昊 on 2023/3/8.
//

#include "LogFile.h"

bool LogFile::testOrCreateDir(const string& pathName){
    size_t found = pathName.find_last_of('/');
    string folder = pathName.substr(0, found + 1);
    ifstream fout (folder);
    if(!fout.good()){
        string commend = "mkdir -p " + folder;
        cout << commend << endl;
        system(commend.c_str());
        fout.open(folder.c_str(), ios_base::in);
        return fout.good();
    }
    return true;
}

LogFile::LogFile(string filePath) {
    logFilePath = filePath;
}

void LogFile::addContent(string &content, bool screenPrint) {
    ofstream outfile;
    outfile.open(logFilePath, ios_base::app | ios_base::out);
    if (outfile.is_open()) {
        outfile << content << endl;
        outfile.close();
    }
    else {
        cout << "[ERROR] Cannot open the log file " << logFilePath << " !\n";
        exit(0);
    }
    if(screenPrint){
        printf("%s\n", content.c_str());
    }
    content.clear();
}

string LogFile::getLogFileName(const string &outputFolder, const string &parametersStr, const string &LDPMethod) {
    string logFilepath = outputFolder + LDPMethod + "/" + "logs/log_" + parametersStr;
    logFilepath.append(".log");
    return logFilepath;
}