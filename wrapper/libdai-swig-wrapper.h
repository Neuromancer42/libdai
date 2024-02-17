//
// Created by Yifan Chen on 2022/6/27.
//

#ifndef LIBDAI_LIBDAI_SWIG_WRAPPER_H
#define LIBDAI_LIBDAI_SWIG_WRAPPER_H

#include <dai/alldai.h>
#include <dai_ext/multithread_emalg.h>
#include <dai_ext/causal_fg.h>
#include <dai_ext/causal_bp.h>
#include <dai_ext/causal_em.h>
#include <chrono>

class LibDAISWIGFactorGraph {
    dai::FactorGraph* fg = nullptr;
    dai::InfAlg* infalg = nullptr;
    dai::PropertySet opts;
    bool activated;
    dai::MultiEMAlg* em = nullptr;
    dai::Evidence* evi = nullptr;
    dai::CausalFactorGraph* causal_fg = nullptr;
    dai::CausalEM* causal_em = nullptr;

public:
    static size_t verbose;
    LibDAISWIGFactorGraph(const std::string& fgFileName, int maxiter, int maxtime, double tol, const std::string& alg) throw(std::runtime_error) {
        std::clog << "LibDAI: Loading factor graph from " << fgFileName << std::endl;
        std::string suffix = fgFileName.substr(fgFileName.find_last_of('.') + 1);
        if (suffix == "causal_fg") {
            causal_fg = new dai::CausalFactorGraph();
            causal_fg->ReadFromFile(fgFileName.c_str());
        } else {
            fg = new dai::FactorGraph();
            fg->ReadFromFile(fgFileName.c_str());
        }
        std::clog << "LibDAI: Loaded factor graph from " << fgFileName << std::endl;
        
        auto startTime = std::chrono::steady_clock::now();

        std::clog << "LibDAI: Initializing " << alg << " ("
                  << "maxiter: " << maxiter << ", "
                  << "maxtime: " << maxtime << "sec" << ", "
                  << "tol: " << tol << ")"
                  << std::endl;

        opts.set("maxiter", static_cast<size_t>(maxiter));
        opts.set("maxtime", dai::Real(maxtime));
        opts.set("tol", dai::Real(tol));
        opts.set("verbose", verbose);
        if (alg == "BP") {
            opts.set("updates", std::string("SEQRND"));
            opts.set("logdomain", true);
            if (fg != nullptr)
                infalg = new dai::BP(*fg, opts);
            else if (causal_fg != nullptr)
                infalg = new dai::CausalBP(*causal_fg, opts);
            else
                DAI_THROW(INVALID_FACTORGRAPH_FILE);
        } else if (alg == "MF") {
            opts.set("init", std::string("UNIFORM"));
            opts.set("updates", std::string("NAIVE"));
            
            if (fg != nullptr)
                infalg = new dai::MF(*fg, opts);
            else
                DAI_THROW(INVALID_FACTORGRAPH_FILE);
        }
        
        infalg->init();
        auto initTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: " << alg << " initialized in "
                  << std::chrono::duration_cast<std::chrono::seconds>(initTime - startTime).count()
                  << "sec" << std::endl;
        em = nullptr;
        evi = nullptr;
        activated = false;
    }

    ~LibDAISWIGFactorGraph() throw(std::runtime_error) {
        delete fg;
        delete causal_fg;
        delete infalg;
        delete em;
        delete causal_em;
        delete evi;
        std::clog << "LibDAI: factor graph released" << std::endl;
    }

    void reset() throw(std::runtime_error) {
        std::clog << "LibDAI: Re-Initializing inference engine." << std::endl;
        infalg->init();
        std::clog << "LibDAI: inference engine re-initialized" << std::endl;

        delete em;
        em = nullptr;
        delete causal_em;
        causal_em = nullptr;
        delete evi;
        evi = nullptr;

        activated = false;
    }

    void infer() throw(std::runtime_error) {
        std::clog << "LibDAI: Inference started" << std::endl;
        auto startTime = std::chrono::steady_clock::now();
        double yetToConvergeFraction = infalg->run();
        auto endTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: Inference finished ("
                  << "iterations: " << infalg->Iterations() << ", "
                  << "converge fraction: " << yetToConvergeFraction << ", "
                  << "duration: " << std::chrono::duration_cast<std::chrono::minutes>(endTime - startTime).count() << "min" << ")"
                  << std::endl;
        activated = true;
    }

    void observeFromFile(const std::string& obsFileName) throw(std::runtime_error) {
        std::clog << "LibDAI: Loading observation from " << obsFileName << std::endl;
        std::ifstream obsFile(obsFileName);
        int varIndex;
        bool varValue;
        int cnt = 0;
        while (obsFile >> varIndex >> varValue) {
            observeBernoulli(varIndex, varValue);
            ++cnt;
        }
        std::clog << "LibDAI: " << cnt << " observations loaded." << std::endl;
    }

    void observeBernoulli(int varIndex, bool varValue) throw(std::runtime_error) {
        if (activated) {
            std::cerr << "LibDAI: inference has been activated before, observation takes no effect" << std::endl;
        }
        infalg->clamp(varIndex, varValue ? 1 : 0);
    }

    std::vector<double> queryPostParamFactor(int paramIndex) throw(std::runtime_error) {
        if (!activated) {
            std::cerr << "LibDAI: inference not activated yet, returning init values." << std::endl;
        }
        dai::Var v;
        if (fg != nullptr) {
            v = fg->var(paramIndex);
        } else if (causal_fg != nullptr) {
            v = causal_fg->var(paramIndex);
        }
        dai::TFactor<dai::Real> factor = infalg->belief(v);
        std::vector<double> weights;
        for (size_t i = 0; i < factor.nrStates(); ++i) {
            weights.push_back(factor.get(i));
        }
        return weights;
    }

    // Bernoulli distribution has only one parameter {P(X=1)=p, P(X=0)=1-p}
    double queryPostBernoulliParam(int paramIndex) throw(std::runtime_error) {
        auto ans = queryPostParamFactor(paramIndex);
        return ans[1];
    }

    double getPriorBernoulliParam(int paramIndex) throw(std::runtime_error) {
        if (fg != nullptr)
            return infalg->fg().factor(paramIndex)[1];
        else if (causal_fg != nullptr) {
            dai::CausalFactor f = dynamic_cast<dai::CausalBP*>(infalg)->factor(paramIndex);
            assert(f.type == dai::CausalFactor::Singleton);
            return f.prob();
        } else {
            return -1;
        }
    }

    void initEM(const std::string& evidenceFile, const std::string& emFile) throw(std::runtime_error) {
        if (activated) {
            std::clog << "LibDAI: reset to run EM" << std::endl;
            activated = false;
        }
        if (em != nullptr || causal_em != nullptr) {
            std::clog << "LibDAI: reset previous EM settings" << std::endl;
            delete em;
            em = nullptr;
            delete causal_em;
            causal_em = nullptr;
        }
        // ref: see example from https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/example_sprinkler_em_8cpp-example.html
        auto startTime = std::chrono::steady_clock::now();
        // 1. construct evidence file
        std::clog << "LibDAI: loading evidence from " << evidenceFile << std::endl;
        std::ifstream iEvis(evidenceFile);
        evi = new dai::Evidence();
        {
            std::map<std::string, dai::Var> varMap;
            if (fg != nullptr) {
                for (const auto &v: fg->vars()) {
                    std::stringstream s;
                    s << v.label();
                    varMap[s.str()] = v;
                }
            } else {
                for (const auto &v: causal_fg->vars()) {
                    std::stringstream s;
                    s << v.label();
                    varMap[s.str()] = v;
                }
            }
            evi->addEvidenceTabFile(iEvis, varMap);
        }
        auto tabLoadTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: evidence of " << evi->nrSamples() << " samples loaded in "
                  << std::chrono::duration_cast<std::chrono::seconds>(tabLoadTime - startTime).count()
                  << "sec" << std::endl;

        // 2. load EM spec
        std::clog << "LibDAI: EM started initialization"<< std::endl;
        std::ifstream iEMs(emFile);
        if (auto * causal_bp = dynamic_cast<dai::CausalBP*>(infalg)) {
            causal_em = new dai::CausalEM(*evi, *causal_bp, iEMs);
        } else {
            em = new dai::MultiEMAlg(*evi, *infalg, iEMs);
        }
        auto emInitTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: EM initialized in "
                  << std::chrono::duration_cast<std::chrono::seconds>(emInitTime - tabLoadTime).count()
                  << "sec" << std::endl;
    }

    bool isEMconverged() throw(std::runtime_error) {
        if (em != nullptr)
            return em->hasSatisfiedTermConditions();
        else if (causal_em != nullptr)
            return causal_em->hasSatisfiedTermConditions();
        else {
            std::clog << "LibDAI: EM has not been initialized before check" << std::endl;
            return false;
        }
    }

    void iterateEM(size_t max_jobs) throw(std::runtime_error) {
        if (em == nullptr && causal_em == nullptr) {
            std::clog << "LibDAI: EM has not been initialized before iterate, do nothing" << std::endl;
            return;
        }
        auto iterStartTime = std::chrono::steady_clock::now();
        dai::Real likelihood;
        size_t iter;
        if (em != nullptr) {
            em->setMaxJobs(max_jobs);
            likelihood = em->iterate();
            iter = em->Iterations();
        } else {
            causal_em->setMaxJobs(max_jobs);
            likelihood = causal_em->iterate();
            iter = causal_em->Iterations();
        }
        activated = true;
        auto iterEndTime = std::chrono::steady_clock::now();
        std::clog << "Iteration " << iter
                  << " likelihood: " << likelihood
                  << " time: " << std::chrono::duration_cast<std::chrono::seconds>(iterEndTime - iterStartTime).count() << "sec"
                  << std::endl;
        
    }

    void runEM(const std::string& evidenceFile, const std::string& emFile, size_t max_jobs) throw(std::runtime_error) {
        initEM(evidenceFile, emFile);
        // iterating EM until convergence
        while (!isEMconverged()) {
            iterateEM(max_jobs);
        }
    }

    void dumpVars() throw(std::runtime_error) {
        std::clog << infalg->beliefs() << std::endl;
    }
};

size_t LibDAISWIGFactorGraph::verbose = 0;

#endif //LIBDAI_LIBDAI_SWIG_WRAPPER_H
