//
// Created by Yifan Chen on 2022/6/27.
//

#include "libdai-swig-wrapper.h"

using namespace std;

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars, const vector<int> &pVars);

int main(int argc, char *argv[]) {
    if (argc != 5) {
        cerr << "Usage: ./wrapper-test-bp <factor-graph-file> <obs-file> <query-file> <param-file>" << endl;
        return 1;
    }
    char *fgFileName = argv[1];
    clog << "Loading factor graph " << fgFileName << endl;
    LibDAISWIGFactorGraph fg(fgFileName, 10000000, 10800, 1e-6);
    clog << "Factor graph loaded." << endl;

    char *oFileName = argv[2];
    ifstream oFile(oFileName);
    vector<int> oVars;
    vector<bool> oVals;
    clog << "Loading observation from " << oFileName << endl;
    int oVar;
    bool oVal;
    while (oFile >> oVar >> oVal) {
        oVars.push_back(oVar);
        oVals.push_back(oVal);
    }

    clog << "Observations:";
    for (int i = 0; i < oVars.size(); ++i) {
        clog << " " << (oVals[i] ? "+" : "-") << oVars[i];
    }
    clog << "." << endl;

    char *qFileName = argv[3];
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

    char *pFileName = argv[4];
    ifstream pFile(pFileName);
    vector<int> pVars;
    clog << "Reading parameter variables from " << pFileName << endl;
    int pVar;
    while (pFile >> pVar) {
        pVars.push_back(pVar);
    }
    clog << "Parameter variables:";
    for (int v : pVars) {
        clog << " " << v;
    }
    clog << "." << endl;

//    clog << "Init: " << endl;
//    dumpQueries(fg, qVars, pVars);
//    clog << endl;

    clog << "Prior: " << endl;
    fg.runBP();
    dumpQueries(fg, qVars, pVars);
    clog << endl;

    if(!oVars.empty()) {
        clog << "Posterior: " << endl;
        fg.resetBP();
        for (int i = 0; i < oVars.size(); ++i) {
            fg.observeBernoulli(oVars[i], oVals[i]);
        }
        fg.runBP();
        dumpQueries(fg, qVars, pVars);
    }
    return 0;
}

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars, const vector<int> &pVars) {
    clog << "Q:";
    for (int i = 0; i < qVars.size(); ++i) {
        clog << "\t<" << qVars[i] << "," << fg.queryBernoulliParam(qVars[i]) << ">";
    }
    clog << endl;
    clog << "P:";
    for (int i = 0; i < pVars.size(); ++i) {
        clog << "\t<" << pVars[i];
        vector<double> weights = fg.queryParamFactor(pVars[i]);
        for (int j = 0; j < weights.size(); ++j) {
            clog << "," << weights[j];
        }
        clog << ">";
    }
    clog << endl;
}

