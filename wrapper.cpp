#include <dai/alldai.h>
#include <dai_ext/causal_bp.h>
#include <dai_ext/causal_fg.h>
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

static bool is_causal;
static FactorGraph fg;
static CausalFactorGraph causal_fg;
static PropertySet opts;
static unique_ptr<BP> bp;
static unique_ptr<CausalBP> causal_bp;
static map<int, bool> clamps;

void initBP() {
    if (is_causal)
        causal_bp.reset(new CausalBP(causal_fg, opts));
    else
        bp.reset(new BP(fg, opts));
    for (const auto& clamp : clamps) {
        int varIndex = clamp.first;
        bool varValue = clamp.second;
        if (is_causal)
            causal_bp->clamp(varIndex, varValue ? 1 : 0);
        else
            bp->clamp(varIndex, varValue ? 1 : 0);
    }
    if (is_causal)
        causal_bp->init();
    else
        bp->init();
}

void queryVariable() {
    int varIndex;
    cin >> varIndex;
    clog << __LOGSTR__ << "Q " << varIndex << endl;

    // auto ans = bp.belief(fg.var(varIndex)).get(1);
    Real ans = is_causal ? causal_bp->newBelief(varIndex)
                    : bp->newBelief(varIndex);
    clog << __LOGSTR__ << "Returning " << ans << "." << endl;
    cout << ans << endl;
}

void queryFactor() {
    int factorIndex, valueIndex;
    cin >> factorIndex >> valueIndex;
    clog << __LOGSTR__ << "FQ " << factorIndex << " " << valueIndex << endl;

    auto ans = is_causal ? causal_bp->beliefF(factorIndex).get(valueIndex)
                    : bp->beliefF(factorIndex).get(valueIndex); 
    clog << __LOGSTR__ << "Returning " << ans << "." << endl;
    cout << ans << endl;
}

void runBP() {
    double tolerance;
    size_t minIters, maxIters, histLength;
    cin >> tolerance >> minIters >> maxIters >> histLength;
    clog << __LOGSTR__ << "BP " << tolerance << " " << minIters << " " << maxIters << " " << histLength << endl;

    initBP();
    double yetToConvergeFraction = 
            is_causal ? causal_bp->run(tolerance, minIters, maxIters, histLength)
                    : bp->run(tolerance, minIters, maxIters, histLength);
    cout << yetToConvergeFraction << endl;
}

void clamp() {
    int varIndex;
    string varValueStr;
    cin >> varIndex >> varValueStr;
    assert(varValueStr == "true" || varValueStr == "false");
    clog << __LOGSTR__ << "O " << varIndex << " " << varValueStr << endl;

    bool varValue = (varValueStr == "true");
    clamps[varIndex] = varValue;
    cout << "O " << varIndex << " " << varValueStr << endl;
}

void unclamp() {
    int varIndex;
    cin >> varIndex;
    assert(clamps.find(varIndex) != clamps.end());
    clamps.erase(varIndex);
    clog << __LOGSTR__ << "UC " << varIndex << endl;
    cout << "UC " << varIndex << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << __LOGSTR__ << "Insufficient number of arguments." << endl;
        return 1;
    }
    string factorGraphFileName = argv[1];
    is_causal = factorGraphFileName.substr(
            factorGraphFileName.find_last_of(".") + 1
            ) == "causal_fg";
    clog << __LOGSTR__ << "Hello!" << endl
                       << "wrapper.cpp compiled on " << __DATE__ << " at " << __TIME__ << "." << endl;
    if (is_causal)
        causal_fg.ReadFromFile(factorGraphFileName.c_str());
    else
        fg.ReadFromFile(factorGraphFileName.c_str());
    clog << __LOGSTR__ << "Finished reading factor graph." << endl;

    opts.set("maxiter", static_cast<size_t>(10000000));
    opts.set("maxtime", Real(57600));
    opts.set("tol", Real(1e-6));
    opts.set("verb", static_cast<size_t>(1));
    if (is_causal)
        opts.set("updates", string("PARALL"));
    else
        opts.set("updates", string("SEQRND")); // "SEQRND", or "PARALL", or "SEQFIX", or "SEQRNDPAR"
    opts.set("logdomain", true);

    initBP();

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
        } else if (cmdType == "UC") {
            unclamp();
        } else {
            assert(cmdType == "NL");
            cout << endl;
        }
    }

    clog << __LOGSTR__ << "Bye!" << endl;
    return 0;
}
