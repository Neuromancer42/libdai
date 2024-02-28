//
// Created by Yifan Chen on 2022/6/27.
//

#include "libdai-swig-wrapper.h"

using namespace std;

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars, const vector<int> &pVars);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: ./wrapper-test-bp <factor-graph-file> <query-file> [<param-file> [<obs-file>]]" << endl;
        return 1;
    }
    LibDAISWIGFactorGraph::verbose = 10;
    char *fgFileName = argv[1];
    clog << "Loading factor graph " << fgFileName << endl;
    LibDAISWIGFactorGraph fg(fgFileName, 10000, 3600, 1e-6, "BP");
    clog << "Factor graph loaded." << endl;

//    clog << "Init: " << endl;
//    dumpQueries(fg, qVars, pVars);
//    clog << endl;

    clog << "Running Prior: " << endl;
    fg.infer();

    if (argc < 3) {
        cerr << "No queries, quit" << endl;
        return 0;
    }
    char *qFileName = argv[2];
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


    vector<int> pVars;
    if (argc < 4) {
        cerr << "No params specified, use empty" << endl;
    } else {
        char *pFileName = argv[3];
        ifstream pFile(pFileName);
        clog << "Reading parameter variables from " << pFileName << endl;
        int pVar;
        while (pFile >> pVar) {
            pVars.push_back(pVar);
        }
        clog << "Parameter variables:";
        for (int v: pVars) {
            clog << " " << v;
        }
        clog << "." << endl;
    }

    clog << "Prior: " << endl;
    dumpQueries(fg, qVars, pVars);
    clog << endl;

    if (argc < 5) {
        cerr << " no obs specified, return" << endl;
        return 0;
    }
    char *oFileName = argv[4];
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
    if(!oVars.empty()) {
        clog << "Posterior: " << endl;
        fg.reset();
        for (int i = 0; i < oVars.size(); ++i) {
            fg.observeBernoulli(oVars[i], oVals[i]);
        }
        fg.infer();
        dumpQueries(fg, qVars, pVars);
    }
    return 0;
}

void dumpQueries(LibDAISWIGFactorGraph &fg, const vector<int> &qVars, const vector<int> &pVars) {
    clog << "Q:";
    for (int i = 0; i < qVars.size(); ++i) {
        clog << "\t<" << qVars[i] << "," << fg.queryPostBernoulliParam(qVars[i]) << ">";
    }
    clog << endl;
    clog << "P:";
    for (int i = 0; i < pVars.size(); ++i) {
        clog << "\t<" << pVars[i];
        vector<double> weights = fg.queryPostParamFactor(pVars[i]);
        for (int j = 0; j < weights.size(); ++j) {
            clog << "," << weights[j];
        }
        clog << ">";
    }
    clog << endl;
}

