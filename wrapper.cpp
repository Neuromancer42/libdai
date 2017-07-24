#include <dai/alldai.h>
using namespace dai;

#include <boost/lexical_cast.hpp>

#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
using namespace std;

static size_t numIters;
static FactorGraph fg;
static map<int, bool> evidence;
static BP *bp = nullptr, *bpRestore = nullptr;
static bool stale = true;

bool runBP() {
    if (stale) {
        clog << __LOGSTR__ << "Computing marginal probabilities." << endl;
        bool converged = bp->run(numIters);
        if (converged) {
            bpRestore = new BP(*bp);
            clog << __LOGSTR__ << "Finished computing marginal probabilities. " << endl
                               << "Executed " << bp->Iterations() << " iterations." << endl;
        } else {
            assert(bpRestore != nullptr);
            bp = new BP(*bpRestore);
            clog << __LOGSTR__ << "Belief propagation did not converge after " << numIters << " iterations. " << endl
                               << "Restoring to backup point." << endl;
        }
        stale = false;
        return converged;
    } else {
        return true;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << __LOGSTR__ << "Insufficient number of arguments." << endl;
        return 1;
    }
    char *factorGraphFileName = argv[1];
    double tolerance = boost::lexical_cast<double>(argv[2]);
    numIters = boost::lexical_cast<size_t>(argv[3]);

    clog << __LOGSTR__ << "Hello!" << endl
                       << "wrapper.cpp compiled on " << __DATE__ << " at " << __TIME__ << "." << endl;
    fg.ReadFromFile(factorGraphFileName);
    clog << __LOGSTR__ << "Finished reading factor graph." << endl;

    PropertySet opts;
    opts.set("maxiter", static_cast<size_t>(10000000));
    opts.set("tol", Real(tolerance));
    opts.set("verb", static_cast<size_t>(1));
    opts.set("updates", string("SEQRND")); // "PARALL", or "SEQFIX"
    opts.set("logdomain", false);
    // opts.set("damping", string("0.2"));

    bp = new BP(fg, opts);
    bp->init();
    bpRestore = new BP(bp);
    stale = true;

    string cmdType;
    while (cin >> cmdType) {
        clog << __LOGSTR__ << "Read command " << cmdType << endl;
        if (cmdType == "Q") {
            int varIndex;
            cin >> varIndex;
            clog << __LOGSTR__ << "Q " << varIndex << endl;
            if (evidence.find(varIndex) == evidence.end()) {
                clog << __LOGSTR__ << "Variable not part of evidence. Querying LibDAI." << endl;
                cout << bp->belief(fg.var(varIndex)).get(1) << endl;
            } else {
                clog << __LOGSTR__ << "Variable previously clamped as evidence. Immediately returning value." << endl;
                cout << evidence.at(varIndex) ? 1 : 0 << endl;
            }
        } else if (cmdType == "FQ") {
            int factorIndex, valueIndex;
            cin >> factorIndex >> valueIndex;
            clog << __LOGSTR__ << "FQ " << factorIndex << " " << valueIndex << endl;
            cout << bp->beliefF(factorIndex).get(valueIndex) << endl;
        } else if (cmdType == "BP") {
            bool converged = runBP();
            if (converged) {
                cout << "converged" << endl;
            } else {
                cout << "diverged" << endl;
            }
        } else if (cmdType == "O") {
            int varIndex;
            string varValueStr;
            cin >> varIndex >> varValueStr;
            assert(varValueStr == "true" || varValueStr == "false");
            bool varValue = (varValueStr == "true");
            clog << __LOGSTR__ << "Clamping variable " << varIndex << " to value " << varValue << "." << endl;
            bp->clamp(varIndex, varValue ? 1 : 0);
            stale = true;
        } else {
            assert(cmdType == "NL");
            cout << endl;
        }
    }

    clog << __LOGSTR__ << "Bye!" << endl;
    delete bp;
    delete bpRestore;
    return 0;
}
