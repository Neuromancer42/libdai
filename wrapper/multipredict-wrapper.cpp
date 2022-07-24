#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
using namespace std;

#include <dai/alldai.h>
using namespace dai;


static FactorGraph fg;
static PropertySet opts;
static BP bp;
static map<int, bool> clamps;

TFactor<Real> queryParam(int paraIndex) {
    clog << "Query factor " << paraIndex << endl;
    auto ans = bp.belief(fg.var(paraIndex));
    clog << "Returning " << ans << "." << endl;
    return ans;
}

void observe(int varIndex, bool varValue) {
    clog << "Observe var " << varIndex << " " << varValue << endl;
    clamps[varIndex] = varValue;
}

void observeFromFile(const string& obsFileName) {
    clog << "Loading observation from " << obsFileName << endl;
    ifstream obsFile(obsFileName);
    int varIndex;
    bool varValue;
    int cnt = 0;
    while(obsFile >> varIndex >> varValue) {
        observe(varIndex, varValue);
        ++cnt;
    }
    clog << cnt << " observations loaded." << endl;
}

void initBP() {
    bp = BP(fg, opts);
    for (const auto& clamp: clamps) {
        int varIndex = clamp.first;
        bool varValue = clamp.second;
        bp.clamp(varIndex, varValue ? 1 : 0);
    }
    bp.init();
}

double runBP(double tolerance, size_t minIters, size_t maxIters, size_t histLength) {
    clog << "BP " << tolerance << " " << minIters << "~" << maxIters << " " << histLength << endl;
    //double yetToConvergeFraction  = bp.run(tolerance, minIters, maxIters, histLength);
    double yetToConvergeFraction = bp.run();
    clog << "BP finished, " << yetToConvergeFraction << endl;
    return yetToConvergeFraction;
}

int main(int argc, char *argv[]) {
    if (argc != 4 && argc != 5 && argc != 6 && argc != 7) {
        cerr << "Usage: ./wrapper <factor-graph> <obs-file> <query-file> [<prediction-file> [<param-file> [<weights-file>]]]" << endl;
        return 1;
    }

    char *fgFileName = argv[1];
    clog << "Loading factor graph " << fgFileName << endl;
    fg.ReadFromFile(fgFileName);
    clog << "Factor graph loaded." << endl;

    char *obsFileName = argv[2];
    observeFromFile(obsFileName);

    opts.set("maxiter", static_cast<size_t>(10000000));
    opts.set("maxtime", Real(3600));
    opts.set("tol", Real(1e-6));
    opts.set("updates", string("SEQRND"));
    opts.set("logdomain", true);

    initBP();

    double tolerance = 1e-6;
    size_t minIter = 500;
    size_t maxIter = 1000;
    size_t histLength = 100;
    runBP(tolerance, minIter, maxIter, histLength);

    char *qFileName = argv[3];
    ifstream qFile(qFileName);
    int varIndex;
    clog << "Querying variables" << endl;
    ostringstream oss;
    while (qFile >> varIndex) {
        const TFactor<Real>& paraFactor = queryParam(varIndex);
        oss << varIndex;
        oss << "\t";
        oss << setprecision(numeric_limits<double>::max_digits10)
            << paraFactor.get(1);
        oss << endl;
    }

    if (argc >= 5) {
        char *pFileName = argv[4];
        clog << "Dumping prediction to " << pFileName << endl;
        ofstream pFile(pFileName);
        pFile << oss.str();
    } else {
        cout << oss.str();
    };

    if (argc >= 6) {
        char *pFileName = argv[5];
        ifstream pFile(pFileName);
        int paraIndex;
        clog << "Querying parameters" << endl;
        ostringstream oss;
        while (pFile >> paraIndex) {
            const TFactor<Real>& paraFactor = queryParam(paraIndex);
            oss << paraIndex;
            for (size_t i = 0; i < paraFactor.nrStates(); ++i) {
                oss << "\t";
                oss << setprecision(numeric_limits<double>::max_digits10)
                    << paraFactor.get(i);
            }
            oss << endl;
        }
        if (argc == 7) {
            char *wFileName = argv[6];
            clog << "Dumping updated parameters to " << wFileName << endl;
            ofstream wFile(wFileName);
            wFile << oss.str();
        } else {
            cout << oss.str();
        };
    }

    clog << "All: " << endl;
    clog << bp.beliefs() << endl;

    return 0;
}