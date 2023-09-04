//
// Created by Yifan Chen on 2022/6/27.
//

#ifndef LIBDAI_LIBDAI_SWIG_WRAPPER_H
#define LIBDAI_LIBDAI_SWIG_WRAPPER_H

#include <dai/alldai.h>
#include <dai_ext/multithread_emalg.h>
#include <chrono>

class LibDAISWIGFactorGraph {
    dai::FactorGraph* fg;
    dai::InfAlg* infalg;
    dai::PropertySet opts;
    bool activated;
    dai::MultiEMAlg* em;
    dai::Evidence* evi;

public:
    LibDAISWIGFactorGraph(const std::string& fgFileName, int maxiter, int maxtime, double tol, const std::string& alg) {
        //std::clog << "LibDAI: Loading factor graph from " << fgFileName << std::endl;
        fg = new dai::FactorGraph();
        fg->ReadFromFile(fgFileName.c_str());
        std::clog << "LibDAI: Loaded factor graph from " << fgFileName << std::endl;

        opts.set("maxiter", static_cast<size_t>(maxiter));
        opts.set("maxtime", dai::Real(maxtime));
        opts.set("tol", dai::Real(tol));

        if (alg == "BP") {
            std::clog << "LibDAI: Initializing loopy BP ("
                      << "maxiter: " << opts.getStringAs<size_t>("maxiter") << ", "
                      << "maxtime: " << opts.getStringAs<dai::Real>("maxtime") << "sec" << ", "
                      << "tol: " << opts.getStringAs<dai::Real>("tol") << ")"
                      << std::endl;
            auto startTime = std::chrono::steady_clock::now();

            opts.set("updates", std::string("SEQFIX"));
            opts.set("logdomain", true);
            infalg = new dai::BP(*fg, opts);
            infalg->init();

            auto initTime = std::chrono::steady_clock::now();
            std::clog << "LibDAI: BP initialized in "
                      << std::chrono::duration_cast<std::chrono::seconds>(initTime - startTime).count()
                      << "sec" << std::endl;
        } else if (alg == "MF") {
            std::clog << "LibDAI: Initializing MF ("
                      << "maxiter: " << opts.getStringAs<size_t>("maxiter") << ", "
                      << "maxtime: " << opts.getStringAs<dai::Real>("maxtime") << "sec" << ", "
                      << "tol: " << opts.getStringAs<dai::Real>("tol") << ")"
                      << std::endl;
            auto startTime = std::chrono::steady_clock::now();

            opts.set("init", std::string("UNIFORM"));
            opts.set("updates", std::string("NAIVE"));
            infalg = new dai::MF(*fg, opts);
            infalg->init();

            auto initTime = std::chrono::steady_clock::now();
            std::clog << "LibDAI: MF initialized in "
                      << std::chrono::duration_cast<std::chrono::seconds>(initTime - startTime).count()
                      << "sec" << std::endl;
        }
        em = nullptr;
        evi = nullptr;
        activated = false;
    }

    ~LibDAISWIGFactorGraph() {
        delete fg;
        delete infalg;
        delete em;
        delete evi;
        std::clog << "LibDAI: factor graph released" << std::endl;
    }

    void reset() {
        std::clog << "LibDAI: Re-Initializing inference engine." << std::endl;
        infalg->init();
        std::clog << "LibDAI: inference engine re-initialized" << std::endl;

        delete em;
        em = nullptr;
        delete evi;
        evi = nullptr;

        activated = false;
    }

    void infer() {
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

    void observeFromFile(const std::string& obsFileName) {
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

    void observeBernoulli(int varIndex, bool varValue) {
        if (activated) {
            std::cerr << "LibDAI: inference has been activated before, observation takes no effect" << std::endl;
        }
        infalg->clamp(varIndex, varValue ? 1 : 0);
    }

    std::vector<double> queryPostParamFactor(int paramIndex) {
        if (!activated) {
            std::cerr << "LibDAI: inference not activated yet, returning init values." << std::endl;
        }
        dai::TFactor<dai::Real> factor = infalg->belief(fg->var(paramIndex));
        std::vector<double> weights;
        for (size_t i = 0; i < factor.nrStates(); ++i) {
            weights.push_back(factor.get(i));
        }
        return weights;
    }

    // Bernoulli distribution has only one parameter {P(X=1)=p, P(X=0)=1-p}
    double queryPostBernoulliParam(int paramIndex) {
        auto ans = queryPostParamFactor(paramIndex);
        return ans[1];
    }

    double getPriorBernoulliParam(int paramIndex) {
        return infalg->fg().factor(paramIndex)[1];
    }

    void initEM(const std::string& evidenceFile, const std::string& emFile) {
        if (activated) {
            std::clog << "LibDAI: reset to run EM" << std::endl;
            activated = false;
        }
        if (em != nullptr) {
            std::clog << "LibDAI: reset previous EM settings" << std::endl;
            delete em;
            em = nullptr;
        }
        // ref: see example from https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/example_sprinkler_em_8cpp-example.html
        auto startTime = std::chrono::steady_clock::now();
        // 1. construct evidence file
        std::clog << "LibDAI: loading evidence from " << evidenceFile << std::endl;
        std::ifstream iEvis(evidenceFile);
        evi = new dai::Evidence();
        evi->addEvidenceTabFile(iEvis, *fg);
        auto tabLoadTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: evidence of " << evi->nrSamples() << " samples loaded in "
                  << std::chrono::duration_cast<std::chrono::seconds>(tabLoadTime - startTime).count()
                  << "sec" << std::endl;

        // 2. load EM spec
        std::clog << "LibDAI: EM started initialization"<< std::endl;
        std::ifstream iEMs(emFile);
        em = new dai::MultiEMAlg(*evi, *infalg, iEMs);
        auto emInitTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: EM initialized in "
                  << std::chrono::duration_cast<std::chrono::seconds>(emInitTime - tabLoadTime).count()
                  << "sec" << std::endl;
    }

    bool isEMconverged() {
        if (em == nullptr) {
            std::clog << "LibDAI: EM has not been initialized before check" << std::endl;
            return false;
        }
        return em->hasSatisfiedTermConditions();
    }

    void iterateEM(size_t max_jobs) {
        if (em == nullptr) {
            std::clog << "LibDAI: EM has not been initialized before iterate, do nothing" << std::endl;
            return;
        }
        auto iterStartTime = std::chrono::steady_clock::now();
        em->setMaxJobs(max_jobs);
        dai::Real I = em->iterate();
        activated = true;
        auto iterEndTime = std::chrono::steady_clock::now();
        std::clog << "Iteration " << em->Iterations()
                  << " likelihood: " << I
                  << " time: " << std::chrono::duration_cast<std::chrono::seconds>(iterEndTime - iterStartTime).count() << "sec"
                  << std::endl;
    }

    void runEM(const std::string& evidenceFile, const std::string& emFile, size_t max_jobs) {
        initEM(evidenceFile, emFile);
        // iterating EM until convergence
        while (!isEMconverged()) {
            iterateEM(max_jobs);
        }
    }
};

#endif //LIBDAI_LIBDAI_SWIG_WRAPPER_H
