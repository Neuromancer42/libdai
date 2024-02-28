//
// Created by Yifan Chen on 2024/1/14.
//
#include <dai_ext/causal_em.h>
#include <future>
#include <thread>

using namespace dai;


void CausalSharedParam::setPermsAndVarSetsFromVarOrders() {
    if( _varorders.size() == 0 )
        return;
    DAI_ASSERT( _estimation != NULL );
    _expectations = new Prob(_estimation->probSize(), 0);

//    // Construct the permutation objects and the varsets
//    for( FactorOrientations::const_iterator foi = _varorders.begin(); foi != _varorders.end(); ++foi ) {
//        VarSet vs;
//        _perms[foi->first] = calculatePermutation( foi->second, vs );
//        _varsets[foi->first] = vs;
//        DAI_ASSERT( (BigInt)_estimation->probSize() == vs.nrStates() );
//    }
}


CausalSharedParam::CausalSharedParam( const FactorOrientations &varorders, ParameterEstimation *estimation, bool ownPE )
        : _varorders(varorders), _estimation(estimation), _ownEstimation(ownPE), _expectations(NULL)
{
    // Calculate the necessary permutations and varsets
    setPermsAndVarSetsFromVarOrders();
}


CausalSharedParam::CausalSharedParam( std::istream &is, const CausalFactorGraph &causal_fg )
        : _varorders(), _estimation(NULL), _ownEstimation(true), _expectations(NULL)
{
    // Read the desired parameter estimation method from the stream
    std::string est_method;
    PropertySet props;
    is >> est_method;
    is >> props;
    DAI_ASSERT(props.getStringAs<size_t>("target_dim") == 2 && props.getStringAs<size_t>("total_dim") == 2);
    
    // Construct a corresponding object
    _estimation = ParameterEstimation::construct( est_method, props );
    // Read in the factors that are to be estimated
    size_t num_factors;
    is >> num_factors;
    for( size_t sp_i = 0; sp_i < num_factors; ++sp_i ) {
        std::string line;
        while( line.size() == 0 && getline(is, line) )
            ;

        std::vector<std::string> fields = tokenizeString( line, true, " \t" );

        // Lookup the factor in the factorgraph
        if( fields.size() < 1 )
            DAI_THROWE(INVALID_EMALG_FILE,"Empty line unexpected");
        std::istringstream iss;
        iss.str( fields[0] );
        size_t factor;
        iss >> factor;
        const CausalFactor & fac = causal_fg.factor(factor);
        DAI_ASSERT(fields.size() == 2 && fac.type == CausalFactor::Singleton);
        _varorders[factor] = fac.head();
        
//        const VarSet & vs = fac.vars();
//        // Construct the vector of Vars
//        std::vector<Var> var_order;
//        var_order.reserve( vs.size() );
//        for( size_t fi = 1; fi < fields.size(); ++fi ) {
//            // Lookup a single variable by label
//            size_t label;
//            std::istringstream labelparse( fields[fi] );
//            labelparse >> label;
//            VarSet::const_iterator vsi = vs.begin();
//            for( ; vsi != vs.end(); ++vsi )
//                if( vsi->label() == label )
//                    break;
//            if( vsi == vs.end() )
//                DAI_THROWE(INVALID_EMALG_FILE,"Specified variables do not match the factor variables");
//            var_order.push_back( *vsi );
//        }
//        _varorders[factor] = var_order;
    }

    // Calculate the necessary permutations
    setPermsAndVarSetsFromVarOrders();
}


void CausalSharedParam::collectExpectations( CausalBP &alg ) {
    for( FactorOrientations::iterator i = _varorders.begin(); i != _varorders.end(); ++i ) {
        Var & v = i->second;
        Factor b = alg.belief(v);
        Prob p( b.nrStates(), 0.0 );
        for( size_t entry = 0; entry < b.nrStates(); ++entry )
            p.set(entry, b[entry]);
//            p.set( entry, b[perm.convertLinearIndex(entry)] ); // apply inverse permutation
        (*_expectations) += p;
    }
}


void CausalSharedParam::setParameters( CausalFactorGraph &causal_fg ) {
    Prob p = _estimation->estimate(this->currentExpectations());
    for( FactorOrientations::iterator i = _varorders.begin(); i != _varorders.end(); ++i ) {
        
        Var& v = i->second;
        
        CausalFactor f( v, p[1]);
        
        causal_fg.setFactor( i->first, f );
    }
}


CausalMaxStep::CausalMaxStep( std::istream &is, const CausalFactorGraph &fg_varlookup ) : _params() {
    size_t num_params = -1;
    is >> num_params;
    _params.reserve( num_params );
    for( size_t i = 0; i < num_params; ++i )
        _params.push_back( CausalSharedParam( is, fg_varlookup ) );
}


void CausalMaxStep::addExpectations( CausalBP &alg ) {
    for( size_t i = 0; i < _params.size(); ++i )
        _params[i].collectExpectations( alg );
}


void CausalMaxStep::maximize( CausalFactorGraph &fg ) {
    for( size_t i = 0; i < _params.size(); ++i )
        _params[i].setParameters( fg );
}


void CausalMaxStep::clear( ) {
    for( size_t i = 0; i < _params.size(); ++i )
        _params[i].clear( );
}

const std::string CausalEM::MAX_ITERS_KEY("max_iters");
const std::string CausalEM::LOG_Z_TOL_KEY("log_z_tol");
const size_t CausalEM::MAX_ITERS_DEFAULT = 30;
const Real CausalEM::LOG_Z_TOL_DEFAULT = 0.01;
const size_t CausalEM::MAX_THREADS=std::thread::hardware_concurrency();

CausalEM::CausalEM( const Evidence &evidence, CausalBP &estep, std::istream &msteps_file )
        : _evidence(evidence), _estep(estep), _msteps(), _iters(0), _lastLogZ(), _max_iters(MAX_ITERS_DEFAULT), _log_z_tol(LOG_Z_TOL_DEFAULT), _max_jobs(MAX_THREADS)
{
    msteps_file.exceptions( std::istream::eofbit | std::istream::failbit | std::istream::badbit );
    size_t num_msteps = -1;
    msteps_file >> num_msteps;
    _msteps.reserve(num_msteps);
    for( size_t i = 0; i < num_msteps; ++i )
        _msteps.push_back( CausalMaxStep( msteps_file, estep ) );
}


void CausalEM::setTermConditions( const PropertySet &p ) {
    if( p.hasKey(MAX_ITERS_KEY) )
        _max_iters = p.getStringAs<size_t>(MAX_ITERS_KEY);
    if( p.hasKey(LOG_Z_TOL_KEY) )
        _log_z_tol = p.getStringAs<Real>(LOG_Z_TOL_KEY);
}


bool CausalEM::hasSatisfiedTermConditions() const {
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


Real CausalEM::iterate( CausalMaxStep &mstep ) {
    Real logZ = 0;
    Real likelihood = 0;

    mstep.clear();

    _estep.run();
    logZ = _estep.logZ();

    // define the task of computing estimations in parallel
    //std::mutex mstep_mutex;
    auto e_func = [](CausalBP * clamped, const Evidence::const_iterator::value_type& e, size_t id) -> Real {
        std::clog << std::string("LibDAI: Expectation task ") + std::to_string(id) + " started\n";
        // apply evidence
        for (const auto & i : e)
            clamped->clamp(clamped->findVar(i.first), i.second);
        clamped->init();
        clamped->run();

        Real single_likelihood = clamped->logZ();
        std::clog << std::string("LibDAI: Expectation task ") + std::to_string(id) + " finished after " + std::to_string(clamped->Iterations()) + " iterations, likelihood " + std::to_string(single_likelihood) + "\n";

        return single_likelihood;
    };

    std::clog << std::string("LibDAI: Running in ") + std::to_string(_max_jobs) + " threads" << std::endl;
    size_t gnum = (_evidence.nrSamples() + _max_jobs - 1) / _max_jobs;

    for (size_t gid = 0; gid < gnum; gid++) {
        std::vector<std::future<Real>> exp_tasks;
        std::vector<std::unique_ptr<CausalBP>> clampeds;
        for (size_t id = gid * _max_jobs; id < (gid + 1) * _max_jobs && id < _evidence.nrSamples(); ++id) {
            std::clog << std::string("LibDAI: Preparing task ") + std::to_string(id) + "\n";
            CausalBP* clamped = _estep.clone();
            clampeds.emplace_back(clamped);
            
            bool group_last = id + 1 == (gid + 1) * _max_jobs || id + 1 == _evidence.nrSamples();
            if (group_last) {
                Real res = e_func(clamped, _evidence.begin()[id], id);
                likelihood += res - logZ;
            } else {
                exp_tasks.push_back(std::async(std::launch::async, e_func, clamped, _evidence.begin()[id], id));
            }
        }
        for (auto &t: exp_tasks)
            likelihood += t.get() - logZ;
        for (auto &clamped: clampeds)
            mstep.addExpectations(*clamped);
    }
    
    // Maximization of parameters
    mstep.maximize( _estep );

    return likelihood;
}

Real CausalEM::iterate() {
    Real likelihood;
    for( size_t i = 0; i < _msteps.size(); ++i )
        likelihood = iterate( _msteps[i] );
    _lastLogZ.push_back( likelihood );
    ++_iters;
    return likelihood;
}


void CausalEM::run() {
    while( !hasSatisfiedTermConditions() )
        iterate();
}