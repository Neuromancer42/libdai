#include <dai/alldai.h>
using namespace dai;

#include <boost/lexical_cast.hpp>

#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;

static FactorGraph fg;
static vector<BP> bpVec;

void queryVariable() {
    int varIndex;
    cin >> varIndex;
    clog << __LOGSTR__ << "Q " << varIndex << endl;

    auto ans = bpVec.back().belief(fg.var(varIndex)).get(1);
    clog << __LOGSTR__ << "Returning " << ans << "." << endl;
    cout << ans << endl;
}

void queryFactor() {
    int factorIndex, valueIndex;
    cin >> factorIndex >> valueIndex;
    clog << __LOGSTR__ << "FQ " << factorIndex << " " << valueIndex << endl;

    auto ans = bpVec.back().beliefF(factorIndex).get(valueIndex);
    clog << __LOGSTR__ << "Returning " << ans << "." << endl;
    cout << ans << endl;
}

void runBP() {
    size_t maxIters;
    cin >> maxIters;
    clog << __LOGSTR__ << "BP " << maxIters << endl;

    bool converged = bpVec.back().run(maxIters);
    if (converged) {
        clog << __LOGSTR__ << "Finished computing marginal probabilities. " << endl
                           << "Executed " << bpVec.back().Iterations() << " iterations so far." << endl;
        cout << "converged" << endl;
    } else {
        clog << __LOGSTR__ << "Belief propagation did not converge after " << maxIters << " iterations." << endl;
        cout << "diverged" << endl;
    }
}

void clamp() {
    int varIndex;
    string varValueStr;
    cin >> varIndex >> varValueStr;
    assert(varValueStr == "true" || varValueStr == "false");
    clog << __LOGSTR__ << "O " << varIndex << " " << varValueStr << endl;

    bool varValue = (varValueStr == "true");
    clog << __LOGSTR__ << "Clamping variable " << varIndex << " to value " << varValue << "." << endl;
    bpVec.back().clamp(varIndex, varValue ? 1 : 0);
    cout << "O " << varIndex << " " << varValueStr << endl;
}

void bpVecPush() {
    clog << __LOGSTR__ << "PUSH" << endl;
    bpVec.push_back(bpVec.back());
    clog << __LOGSTR__ << "New bpVec.size(): " << bpVec.size() << endl;
    cout << bpVec.size() << endl;
}

void bpVecPop() {
    clog << __LOGSTR__ << "POP" << endl;
    bpVec.pop_back();
    assert(bpVec.size() > 0);
    clog << __LOGSTR__ << "New bpVec.size(): " << bpVec.size() << endl;
    cout << bpVec.size() << endl;
}

void bpVecFlatten() {
    clog << __LOGSTR__ << "FLATTEN" << endl;
    bpVec.erase(bpVec.begin(), bpVec.begin() + (bpVec.size() - 1));
    assert(bpVec.size() == 1);
    clog << __LOGSTR__ << "New bpVec.size(): " << bpVec.size() << endl;
    cout << bpVec.size() << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << __LOGSTR__ << "Insufficient number of arguments." << endl;
        return 1;
    }
    char *factorGraphFileName = argv[1];
    double tolerance = boost::lexical_cast<double>(argv[2]);

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

    bpVec.push_back(BP(fg, opts));

    string cmdType;
    while (cin >> cmdType) {
        clog << __LOGSTR__ << "Read command " << cmdType << endl;
        if (cmdType == "Q") {
            queryVariable();
        } else if (cmdType == "FQ") {
            queryFactor();
        } else if (cmdType == "BP") {
            runBP();
        } else if (cmdType == "O") {
            clamp();
        } else if (cmdType == "PUSH") {
            bpVecPush();
        } else if (cmdType == "POP") {
            bpVecPop();
        } else if (cmdType == "FLATTEN") {
            bpVecFlatten();
        } else {
            assert(cmdType == "NL");
            cout << endl;
        }
    }

    clog << __LOGSTR__ << "Bye!" << endl;
    return 0;
}
