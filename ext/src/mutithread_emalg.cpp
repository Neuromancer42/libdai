//
// Created by Yifan Chen on 2023/8/21.
//

#include "dai_ext/multithread_emalg.h"
#include <future>
#include <thread>

namespace dai {

    const std::string MultiEMAlg::MAX_ITERS_KEY("max_iters");
    const std::string MultiEMAlg::LOG_Z_TOL_KEY("log_z_tol");
    const size_t MultiEMAlg::MAX_ITERS_DEFAULT = 30;
    const Real MultiEMAlg::LOG_Z_TOL_DEFAULT = 0.01;
    const size_t MultiEMAlg::MAX_THREADS=std::thread::hardware_concurrency();

    MultiEMAlg::MultiEMAlg( const Evidence &evidence, InfAlg &estep, std::istream &msteps_file)
            : _evidence(evidence), _estep(estep), _msteps(), _iters(0), _lastLogZ(), _max_iters(MAX_ITERS_DEFAULT), _log_z_tol(LOG_Z_TOL_DEFAULT), _max_jobs(MAX_THREADS)
    {
        msteps_file.exceptions( std::istream::eofbit | std::istream::failbit | std::istream::badbit );
        size_t num_msteps = -1;
        msteps_file >> num_msteps;
        _msteps.reserve(num_msteps);
        for( size_t i = 0; i < num_msteps; ++i )
            _msteps.push_back( MaximizationStep( msteps_file, estep.fg() ) );
    }


    void MultiEMAlg::setTermConditions( const PropertySet &p ) {
        if( p.hasKey(MAX_ITERS_KEY) )
            _max_iters = p.getStringAs<size_t>(MAX_ITERS_KEY);
        if( p.hasKey(LOG_Z_TOL_KEY) )
            _log_z_tol = p.getStringAs<Real>(LOG_Z_TOL_KEY);
    }


    bool MultiEMAlg::hasSatisfiedTermConditions() const {
        if( _iters >= _max_iters )
            return true;
        else if( _lastLogZ.size() < 3 )
            // need at least 2 to calculate ratio
            // Also, throw away first iteration, as the parameters may not
            // have been normalized according to the estimation method
            return false;
        else {
            Real current = _lastLogZ[_lastLogZ.size() - 1];
            Real previous = _lastLogZ[_lastLogZ.size() - 2];
            if( previous == 0 )
                return false;
            Real diff = current - previous;
            if( diff < 0 ) {
                std::cerr << "Error: in EM log-likehood decreased from " << previous << " to " << current << std::endl;
                return true;
            }
            return (diff / fabs(previous)) <= _log_z_tol;
        }
    }


    Real MultiEMAlg::iterate( MaximizationStep &mstep ) {
        Real logZ = 0;
        Real likelihood = 0;

        mstep.clear();

        _estep.run();
        logZ = _estep.logZ();

        // define the task of computing estimations in parallel
        //std::mutex mstep_mutex;
        auto e_func = [logZ](InfAlg * clamped, const Evidence::const_iterator::value_type& e, size_t id) -> Real {
            std::clog << std::string("LibDAI: Expectation task ") + std::to_string(id) + " started\n";
            // apply evidence
            for (const auto & i : e)
                clamped->clamp(clamped->fg().findVar(i.first), i.second);
            clamped->init();
            clamped->run();

            Real single_likelihood = clamped->logZ() - logZ;
            std::clog << std::string("LibDAI: Expectation task ") + std::to_string(id) + " finished, likelihood " + std::to_string(single_likelihood) + "\n";

            return single_likelihood;
        };


        std::clog << std::string("LibDAI: Running in ") + std::to_string(_max_jobs) + " threads" << std::endl;
        size_t gnum = (_evidence.nrSamples() + _max_jobs - 1) / _max_jobs;

        for (size_t gid = 0; gid < gnum; gid++) {
            std::vector<std::future<Real>> exp_tasks;
            std::vector<std::unique_ptr<InfAlg>> clampeds;
            for (size_t id = gid * _max_jobs; id < (gid + 1) * _max_jobs && id < _evidence.nrSamples(); ++id) {
                std::clog << std::string("LibDAI: Prepareing task ") + std::to_string(id) + "\n";
                InfAlg* clamped = _estep.clone();
                clampeds.emplace_back(clamped);
                exp_tasks.push_back(std::async(std::launch::async, e_func, clamped, _evidence.begin()[id], id));
            }
            for (auto &t: exp_tasks)
                likelihood += t.get();
            for (auto &clamped: clampeds)
                mstep.addExpectations(*clamped);
        }

        // Maximization of parameters
        mstep.maximize( _estep.fg() );

        return likelihood;
    }


    Real MultiEMAlg::iterate() {
        Real likelihood;
        for( size_t i = 0; i < _msteps.size(); ++i )
            likelihood = iterate( _msteps[i] );
        _lastLogZ.push_back( likelihood );
        ++_iters;
        return likelihood;
    }


    void MultiEMAlg::run() {
        while( !hasSatisfiedTermConditions() )
            iterate();
    }

}