//
// Created by Yifan Chen on 2022/6/27.
//

#ifndef LIBDAI_LIBDAI_SWIG_WRAPPER_H
#define LIBDAI_LIBDAI_SWIG_WRAPPER_H

#include <dai/alldai.h>

class LibDAISWIGFactorGraph {
    dai::FactorGraph* fg;
    dai::BP* bp;
    dai::PropertySet opts;
    bool activated;

public:
    explicit LibDAISWIGFactorGraph(const std::string& fgFileName) {
        //std::clog << "LibDAI: Loading factor graph from " << fgFileName << std::endl;
        fg = new dai::FactorGraph();
        fg->ReadFromFile(fgFileName.c_str());

        //std::clog << "LibDAI: Initializing BP inference engine." << std::endl;
        opts.set("maxiter", static_cast<size_t>(10000000));
        opts.set("maxtime", dai::Real(3600));
        opts.set("tol", dai::Real(1e-6));
        opts.set("updates", std::string("SEQRND"));
        opts.set("logdomain", true);
        bp = new dai::BP(*fg, opts);

        activated = false;
    }

    ~LibDAISWIGFactorGraph() {
        delete fg;
        delete bp;
    }

    void resetBP() {
        //std::clog << "LibDAI: Re-Initializing BP inference engine." << std::endl;
        delete bp;
        bp = new dai::BP(*fg, opts);
        activated = false;
    }

    void runBP() {
        std::clog << "LibDAI: BP started" << std::endl;
        bp->init();
        double yetToConvergeFraction = bp->run();
        std::clog << "LibDAI: BP finished "
            << "(converge fraction: " << yetToConvergeFraction << ")" << std::endl;
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
