//
// Created by Yifan Chen on 2023/8/15.
//

#include "libdai-swig-wrapper.h"

using namespace std;

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars);

int main(int argc, char *argv[]) {
    if (argc != 5) {
        cerr << "Usage: ./wrapper-test-em <factor-graph-file> <evidence-file> <param-file> <query-file>" << endl;
        return 1;
    }
    char *fgFileName = argv[1];
    clog << "Loading factor graph " << fgFileName << endl;
    LibDAISWIGFactorGraph fg(fgFileName, 10000000, 10800, 1e-6);
    clog << "Factor graph loaded." << endl;

    char *eFileName = argv[2];
    clog << "Evidence file: " << eFileName << endl;
    char *emspecFileName = argv[3];
    clog << "EM Spec file: " << emspecFileName << endl;

    char *qFileName = argv[4];
    ifstream qFile(qFileName);
    vector<int> qVars;
    clog << "Reading querying variables from " << qFileName << endl;
    int qVar;
    while (qFile >> qVar) {
        qVars.push_back(qVar);
    }
    clog << "Querying variables:";
    for (int v : qVars) {
        clog << " " << v;
    }
    clog << "." << endl;

//    clog << "Init: " << endl;
//    dumpQueries(fg, qVars, pVars);
//    clog << endl;

    clog << "Prior: " << endl;
    fg.runBP();
    dumpQueries(fg, qVars);
    clog << endl;

    clog << "Posterior: " << endl;
    fg.runEM(eFileName, emspecFileName);
    dumpQueries(fg, qVars);
    return 0;
}

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars) {
    clog << "Q:";
    for (int i = 0; i < qVars.size(); ++i) {
        clog << "\t<" << qVars[i] << "," << fg.queryBernoulliParam(qVars[i]) << ">";
    }
    clog << endl;
}