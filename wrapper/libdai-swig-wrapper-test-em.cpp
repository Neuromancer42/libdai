//
// Created by Yifan Chen on 2023/8/15.
//

#include "libdai-swig-wrapper.h"

using namespace std;

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars);

void dumpParams(LibDAISWIGFactorGraph &fg, const vector<int> &pVars);

int main(int argc, char *argv[]) {
    if (argc != 6) {
        cerr << "Usage: ./wrapper-test-em <factor-graph-file> <evidence-file> <em-file> <query-file> <param-file>" << endl;
        return 1;
    }
    LibDAISWIGFactorGraph::verbose = 10;
    char *fgFileName = argv[1];
    clog << "Loading factor graph " << fgFileName << endl;
    LibDAISWIGFactorGraph fg(fgFileName, 10, 10800, 1e-6, "BP");
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

    char *pFileName = argv[5];
    ifstream pFile(pFileName);
    vector<int> pVars;
    clog << "Reading param variables from " << pFileName << endl;
    int pVar;
    while (pFile >> pVar) {
        pVars.push_back(pVar);
    }
    clog << "Param variables:";
    for (int v : pVars) {
        clog << " " << v;
    }
    clog << "." << endl;

//    clog << "Init: " << endl;
//    dumpQueries(fg, qVars, pVars);
//    clog << endl;

//    dumpParams(fg, pVars);
//    clog << "Prior: " << endl;
//    fg.infer();
//    dumpQueries(fg, qVars);
//    clog << endl;

    clog << "Posterior: " << endl;
    fg.initEM(eFileName, emspecFileName);
    while (!fg.isEMconverged()) {
        fg.iterateEM(1);
        dumpParams(fg, pVars);
    }
    dumpQueries(fg, qVars);
    fg.dumpVars();
    return 0;
}

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars) {
    clog << "Q:";
    for (int qVar : qVars) {
        clog << "\t<" << qVar << "," << fg.queryPostBernoulliParam(qVar) << ">";
    }
    clog << endl;
}

void dumpParams(LibDAISWIGFactorGraph &fg, const vector<int> &pVars) {
    clog << "P:";
    for (int pVar : pVars) {
        clog << "\t<" << pVar << "," << fg.getPriorBernoulliParam(pVar) << ">";
    }
    clog << endl;
}