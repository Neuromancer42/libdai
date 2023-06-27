//
// Created by Yifan Chen on 2022/6/27.
//

#ifndef LIBDAI_LIBDAI_SWIG_WRAPPER_H
#define LIBDAI_LIBDAI_SWIG_WRAPPER_H

#include <dai/alldai.h>
#include <chrono>

class LibDAISWIGFactorGraph {
    dai::FactorGraph* fg;
    dai::BP* bp;
    dai::PropertySet opts;
    bool activated;

public:
    explicit LibDAISWIGFactorGraph(const std::string& fgFileName, int maxiter, int maxtime, double tol) {
        //std::clog << "LibDAI: Loading factor graph from " << fgFileName << std::endl;
        fg = new dai::FactorGraph();
        fg->ReadFromFile(fgFileName.c_str());
        std::clog << "LibDAI: Loaded factor graph from " << fgFileName << std::endl;

        //std::clog << "LibDAI: Initializing BP inference engine." << std::endl;
        opts.set("maxiter", static_cast<size_t>(maxiter));
        opts.set("maxtime", dai::Real(maxtime));
        opts.set("tol", dai::Real(tol));
        opts.set("updates", std::string("SEQRND"));
        opts.set("logdomain", true);
        bp = new dai::BP(*fg, opts);
        activated = false;
        std::clog << "LibDAI: BP initialized ("
            << "maxiter: " << opts.getStringAs<size_t>("maxiter") << ", "
            << "maxtime: " << opts.getStringAs<dai::Real>("maxtime") << "sec" << ", "
            << "tol: " << opts.getStringAs<dai::Real>("tol") << ")"
            << std::endl;
    }

    ~LibDAISWIGFactorGraph() {
        delete fg;
        delete bp;
        std::clog << "LibDAI: factor graph released" << std::endl;
    }

    void resetBP() {
        delete bp;
        std::clog << "LibDAI: Re-Initializing BP inference engine." << std::endl;
        bp = new dai::BP(*fg, opts);
        std::clog << "LibDAI: BP re-initialized" << std::endl;
        activated = false;
    }

    void runBP() {
        std::clog << "LibDAI: BP started" << std::endl;
        auto startTime = std::chrono::steady_clock::now();
        bp->init();
        auto initTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: BP initialized in "
            << std::chrono::duration_cast<std::chrono::seconds>(initTime - startTime).count()
            << "sec" << std::endl;
        double yetToConvergeFraction = bp->run();
        auto endTime = std::chrono::steady_clock::now();
        std::clog << "LibDAI: BP finished ("
            << "iterations: " << bp->Iterations() << ", "
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
        bp->clamp(varIndex, varValue ? 1 : 0);
    }

    std::vector<double> queryParamFactor(int paramIndex) {
        if (!activated) {
            std::cerr << "LibDAI: inference not activated yet, returning init values." << std::endl;
        }
        dai::TFactor<dai::Real> factor = bp->belief(fg->var(paramIndex));
        std::vector<double> weights;
        for (size_t i = 0; i < factor.nrStates(); ++i) {
            weights.push_back(factor.get(i));
        }
        return weights;
    }

    // Bernoulli distribution has only one parameter {P(X=1)=p, P(X=0)=1-p}
    double queryBernoulliParam(int paramIndex) {
        auto ans = queryParamFactor(paramIndex);
        return ans[1];
    }
};

#endif //LIBDAI_LIBDAI_SWIG_WRAPPER_H
