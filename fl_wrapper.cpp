#include <dai/alldai.h>
#include <dai_ext/causal_bp.h>
#include <dai_ext/causal_fg.h>
#include <map>
#include <vector>

using namespace dai;
using namespace std;

dai::FactorGraph* fg = nullptr;
dai::CausalFactorGraph* causal_fg = nullptr;
dai::InfAlg* infalg = nullptr;
dai::PropertySet opts;
#ifdef NDEBUG
size_t verbose = 0;
#else
size_t verbose = 10;
#endif

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "Usage: ./wrapper-test-bp <factor-graph> <ORIG|CAUSAL> <query> [<obs> [<maxiter> [<maxtime> [<tol> [mapping]]]]]" << endl;
        return 1;
    }
    string fgFileName = argv[1];
    std::string suffix = fgFileName.substr(fgFileName.find_last_of('.') + 1);
//    DAI_ASSERT(suffix == "causal_fg");
    
    string type = argv[2];
    DAI_ASSERT(type == "ORIG" || type == "CAUSAL");
    bool orig = type == "ORIG";
    
    vector<int> qVars;
    {
        string qFileName = argv[3];
        ifstream qFile(qFileName);
        clog << "Reading querying variables from " << qFileName << endl;
        int qVar;
        while (qFile >> qVar) {
            qVars.push_back(qVar);
        }
    }
    
    map<int, bool> observation;

    if (argc > 4) {
        char *oFileName = argv[4];
        ifstream oFile(oFileName);
        clog << "Loading observation from " << oFileName << endl;
        int oVar;
        bool oVal;
        while (oFile >> oVar >> oVal) {
            observation.emplace(oVar, oVal);
        }
    }
    
    int maxiter = 10000;
    int maxtime = 3600;
    double tol = 1e-6;
    if (argc > 5) 
        maxiter = std::stoi(argv[5]);
    if (argc > 6)
        maxtime = std::stoi(argv[6]);
    if (argc > 7)
        tol = std::stod(argv[7]);
    
    std::clog << "LibDAI: Loading factor graph from " << fgFileName << std::endl;
    causal_fg = new dai::CausalFactorGraph();
    causal_fg->ReadFromFile(fgFileName.c_str());
    if (orig) {
        std::clog << "LibDAI: Converting to original bp" << std::endl;
        fg = new FactorGraph(causal_fg->factorgraph());
    }
    std::clog << "LibDAI: Loaded factor graph from " << fgFileName << std::endl;

    {
        auto startTime = std::chrono::steady_clock::now();

        std::clog << "LibDAI: Initializing " << type << " ("
                  << "maxiter: " << maxiter << ", "
                  << "maxtime: " << maxtime << "sec" << ", "
                  << "tol: " << tol << ")"
                  << std::endl;

        opts.set("maxiter", static_cast<size_t>(maxiter));
        opts.set("maxtime", dai::Real(maxtime));
        opts.set("tol", dai::Real(tol));
        opts.set("verbose", verbose);
        opts.set("updates", std::string("SEQRND"));
        opts.set("logdomain", true);

        if (orig)
            infalg = new dai::BP(*fg, opts);
        else
            infalg = new dai::CausalBP(*causal_fg, opts);
        infalg->init();

        auto initTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: BP initialized in "
                  << std::chrono::duration_cast<std::chrono::seconds>(initTime - startTime).count()
                  << "sec" << std::endl;
    }

    {
        clog << "Observations:";
        for (auto& obs: observation){
            int oVar = obs.first;
            bool oVal = obs.second;
            clog << " " << (oVal ? "+" : "-") << oVar;
            Var v(oVar, 2);
            size_t varIndex;
            try {
                if (orig) {
                    varIndex = fg->findVar(v);
                } else {
                    varIndex = causal_fg->findVar(v);
                }
            }
            catch (Exception &err) {
                assert(err.code() == Exception::OBJECT_NOT_FOUND);
                clog << "(?)";
                continue;
            }
            infalg->clamp(varIndex, oVal ? 1 : 0);
        }
        clog << "." << endl;
    }
    {
        std::clog << "LibDAI: Inference started" << std::endl;
        auto startTime = std::chrono::steady_clock::now();
        double yetToConvergeFraction = infalg->run();
        auto endTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: Inference finished ("
                  << "iterations: " << infalg->Iterations() << ", "
                  << "converge fraction: " << yetToConvergeFraction << ", "
                  << "duration: " << std::chrono::duration_cast<std::chrono::minutes>(endTime - startTime).count()
                  << "min" << ")"
                  << std::endl;
    }
    map<size_t, Real> probs;
    {
        clog << "Q:";
        for (int qVar : qVars) {
            clog << "\t<" << qVar;
            
            dai::Var v(qVar, 2);
            Real true_belief;
            try {
                true_belief = infalg->belief(v)[1];
            } catch (Exception &err) {
                assert(err.code() == Exception::OBJECT_NOT_FOUND);
                clog << "(?)";
                if (observation.find(qVar) != observation.end()) {
                    true_belief = static_cast<Real>(observation[qVar]);
                } else {
                    true_belief = 0.5;
                }
            }
            clog << "," << true_belief << ">";
            probs.emplace(qVar, true_belief);
        }
        clog << endl;
    }
    
    std::clog << "LibDAI: factor graph released" << std::endl;

    if (argc > 8) {
        string mapFileName = argv[8];
        ifstream mapFile(mapFileName);
        clog << "Loading map from " << mapFileName << endl;
        map<size_t, string> vMap;
        size_t vId;
        string vName;
        while (mapFile >> vId >> vName) {
            vMap.emplace(vId, vName);
        }
        clog << "Map loaded, generating report. " << endl;

        vector<pair<Real, string>> report;
        for (auto & res : probs) {
            size_t qVar = res.first;
            Real true_belief = res.second;
            if (vMap.find(qVar) != vMap.end()) {
                report.emplace_back(true_belief, vMap[qVar]);
            }
        }
        std::sort(report.begin(), report.end());
        for (auto & r : report)
            cout << r.second << " prob_bp = " << r.first << endl;

    } else {
        for (auto &res: probs) {
            size_t qVar = res.first;
            Real true_belief = res.second;
            cout << qVar << "\t" << setw(cout.precision() + 4) << true_belief << endl;
        }
    }
    return 0;
}