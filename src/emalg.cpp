/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2009  Charles Vaske  [cvaske at soe dot ucsc dot edu]
 *  Copyright (C) 2009  University of California, Santa Cruz
 */


#include <dai/util.h>
#include <dai/emalg.h>


namespace dai {


std::map<std::string, ParameterEstimation::ParamEstFactory> *ParameterEstimation::_registry = NULL;


void ParameterEstimation::loadDefaultRegistry() {
    _registry = new std::map<std::string, ParamEstFactory>();
    (*_registry)["ConditionalProbEstimation"] = CondProbEstimation::factory;
}


ParameterEstimation* ParameterEstimation::construct( const std::string &method, const PropertySet &p ) {
    if( _registry == NULL )
        loadDefaultRegistry();
    std::map<std::string, ParamEstFactory>::iterator i = _registry->find(method);
    if( i == _registry->end() )
        DAI_THROWE(UNKNOWN_PARAMETER_ESTIMATION_METHOD, "Unknown parameter estimation method: " + method);
    ParamEstFactory factory = i->second;
    return factory(p);
}


ParameterEstimation* CondProbEstimation::factory( const PropertySet &p ) {
    size_t target_dimension = p.getStringAs<size_t>("target_dim");
    size_t total_dimension = p.getStringAs<size_t>("total_dim");
    Real pseudo_count = 1;
    if( p.hasKey("pseudo_count") )
        pseudo_count = p.getStringAs<Real>("pseudo_count");
    return new CondProbEstimation( target_dimension, Prob( total_dimension, pseudo_count ) );
}


CondProbEstimation::CondProbEstimation( size_t target_dimension, const Prob &pseudocounts )
  : _target_dim(target_dimension), _stats(pseudocounts), _initial_stats(pseudocounts)
{
    DAI_ASSERT( !(_stats.size() % _target_dim) );
}


void CondProbEstimation::addSufficientStatistics( const Prob &p ) {
    _stats += p;
}


Prob CondProbEstimation::estimate() {
    // normalize pseudocounts
    for( size_t parent = 0; parent < _stats.size(); parent += _target_dim ) {
        // calculate norm
        Real norm = 0.0;
        size_t top = parent + _target_dim;
        for( size_t i = parent; i < top; ++i )
            norm += _stats[i];
        if( norm != 0.0 )
            norm = 1.0 / norm;
        // normalize
        for( size_t i = parent; i < top; ++i )
            _stats[i] *= norm;
    }
    // reset _stats to _initial_stats
    Prob result = _stats;
    _stats = _initial_stats;
    return result;
}


Permute SharedParameters::calculatePermutation( const std::vector<Var> &varorder, VarSet &outVS ) {
    // Collect all labels and dimensions, and order them in vs
    std::vector<size_t> dims;
    dims.reserve( varorder.size() );
    std::vector<long> labels;
    labels.reserve( varorder.size() );
    for( size_t i = 0; i < varorder.size(); i++ ) {
        dims.push_back( varorder[i].states() );
        labels.push_back( varorder[i].label() );
        outVS |= varorder[i];
    }

    // Construct the sigma array for the permutation object
    std::vector<size_t> sigma;
    sigma.reserve( dims.size() );
    for( VarSet::iterator set_iterator = outVS.begin(); sigma.size() < dims.size(); ++set_iterator )
        sigma.push_back( find(labels.begin(), labels.end(), set_iterator->label()) - labels.begin() );

    return Permute( dims, sigma );
}


void SharedParameters::setPermsAndVarSetsFromVarOrders() {
    if( _varorders.size() == 0 )
        return;
    DAI_ASSERT( _estimation != NULL );

    // Construct the permutation objects and the varsets
    for( FactorOrientations::const_iterator foi = _varorders.begin(); foi != _varorders.end(); ++foi ) {
        VarSet vs;
        _perms[foi->first] = calculatePermutation( foi->second, vs );
        _varsets[foi->first] = vs;
        DAI_ASSERT( _estimation->probSize() == vs.nrStates() );
    }
}


SharedParameters::SharedParameters( std::istream &is, const FactorGraph &fg_varlookup )
  : _varsets(), _perms(), _varorders(), _estimation(NULL), _deleteEstimation(true)
{
    // Read the desired parameter estimation method from the stream
    std::string est_method;
    PropertySet props;
    is >> est_method;
    is >> props;

    // Construct a corresponding object
    _estimation = ParameterEstimation::construct( est_method, props );

    // Read in the factors that are to be estimated
    size_t num_factors;
    is >> num_factors;
    for( size_t sp_i = 0; sp_i < num_factors; ++sp_i ) {
        std::string line;
        while( line.size() == 0 && getline(is, line) )
            ;

        std::vector<std::string> fields;
        tokenizeString(line, fields, " \t");

        // Lookup the factor in the factorgraph
        if( fields.size() < 1 )
            DAI_THROW(INVALID_EMALG_FILE);
        std::istringstream iss;
        iss.str( fields[0] );
        size_t factor;
        iss >> factor;
        const VarSet &vs = fg_varlookup.factor(factor).vars();
        if( fields.size() != vs.size() + 1 )
            DAI_THROW(INVALID_EMALG_FILE);

        // Construct the vector of Vars
        std::vector<Var> var_order;
        var_order.reserve( vs.size() );
        for( size_t fi = 1; fi < fields.size(); ++fi ) {
            // Lookup a single variable by label
            long label;
            std::istringstream labelparse( fields[fi] );
            labelparse >> label;
            VarSet::const_iterator vsi = vs.begin();
            for( ; vsi != vs.end(); ++vsi )
                if( vsi->label() == label )
                    break;
            if( vsi == vs.end() )
                DAI_THROW(INVALID_EMALG_FILE);
            var_order.push_back( *vsi );
        }
        _varorders[factor] = var_order;
    }

    // Calculate the necessary permutations
    setPermsAndVarSetsFromVarOrders();
}


SharedParameters::SharedParameters( const SharedParameters &sp )
  : _varsets(sp._varsets), _perms(sp._perms), _varorders(sp._varorders), _estimation(sp._estimation), _deleteEstimation(sp._deleteEstimation)
{
    // If sp owns its _estimation object, we should clone it instead
    if( _deleteEstimation )
        _estimation = _estimation->clone();
}


SharedParameters::SharedParameters( const FactorOrientations &varorders, ParameterEstimation *estimation, bool deletePE )
  : _varsets(), _perms(), _varorders(varorders), _estimation(estimation), _deleteEstimation(deletePE)
{
    // Calculate the necessary permutations
    setPermsAndVarSetsFromVarOrders();
}


void SharedParameters::collectSufficientStatistics( InfAlg &alg ) {
    for( std::map< FactorIndex, Permute >::iterator i = _perms.begin(); i != _perms.end(); ++i ) {
        Permute &perm = i->second;
        VarSet &vs = _varsets[i->first];

        Factor b = alg.belief(vs);
        Prob p( b.states(), 0.0 );
        for( size_t entry = 0; entry < b.states(); ++entry )
            p[entry] = b[perm.convert_linear_index(entry)];
        _estimation->addSufficientStatistics( p );
    }
}


void SharedParameters::setParameters( FactorGraph &fg ) {
    Prob p = _estimation->estimate();
    for( std::map<FactorIndex, Permute>::iterator i = _perms.begin(); i != _perms.end(); ++i ) {
        Permute &perm = i->second;
        VarSet &vs = _varsets[i->first];

        Factor f( vs, 0.0 );
        for( size_t entry = 0; entry < f.states(); ++entry )
            f[perm.convert_linear_index(entry)] = p[entry];

        fg.setFactor( i->first, f );
    }
}


void SharedParameters::collectParameters( const FactorGraph &fg, std::vector<Real> &outVals, std::vector<Var> &outVarOrder ) {
    FactorOrientations::iterator it = _varorders.begin();
    if( it == _varorders.end() )
        return;
    FactorIndex I = it->first;
    for( std::vector<Var>::const_iterator var_it = _varorders[I].begin(); var_it != _varorders[I].end(); ++var_it )
        outVarOrder.push_back( *var_it );

    const Factor &f = fg.factor(I);
    DAI_ASSERT( f.vars() == _varsets[I] );
    const Permute &perm = _perms[I];
    for( size_t val_index = 0; val_index < f.states(); ++val_index )
        outVals.push_back( f[perm.convert_linear_index(val_index)] );
}


MaximizationStep::MaximizationStep( std::istream &is, const FactorGraph &fg_varlookup ) : _params() {
    size_t num_params = -1;
    is >> num_params;
    _params.reserve( num_params );
    for( size_t i = 0; i < num_params; ++i )
        _params.push_back( SharedParameters( is, fg_varlookup ) );
}


void MaximizationStep::addExpectations( InfAlg &alg ) {
    for( size_t i = 0; i < _params.size(); ++i )
        _params[i].collectSufficientStatistics( alg );
}


void MaximizationStep::maximize( FactorGraph &fg ) {
    for( size_t i = 0; i < _params.size(); ++i )
        _params[i].setParameters( fg );
}


const std::string EMAlg::MAX_ITERS_KEY("max_iters");
const std::string EMAlg::LOG_Z_TOL_KEY("log_z_tol");
const size_t EMAlg::MAX_ITERS_DEFAULT = 30;
const Real EMAlg::LOG_Z_TOL_DEFAULT = 0.01;


EMAlg::EMAlg( const Evidence &evidence, InfAlg &estep, std::istream &msteps_file )
  : _evidence(evidence), _estep(estep), _msteps(), _iters(0), _lastLogZ(), _max_iters(MAX_ITERS_DEFAULT), _log_z_tol(LOG_Z_TOL_DEFAULT)
{
    msteps_file.exceptions( std::istream::eofbit | std::istream::failbit | std::istream::badbit );
    size_t num_msteps = -1;
    msteps_file >> num_msteps;
    _msteps.reserve(num_msteps);
    for( size_t i = 0; i < num_msteps; ++i )
        _msteps.push_back( MaximizationStep( msteps_file, estep.fg() ) );
}


void EMAlg::setTermConditions( const PropertySet &p ) {
    if( p.hasKey(MAX_ITERS_KEY) )
        _max_iters = p.getStringAs<size_t>(MAX_ITERS_KEY);
    if( p.hasKey(LOG_Z_TOL_KEY) )
        _log_z_tol = p.getStringAs<Real>(LOG_Z_TOL_KEY);
}


bool EMAlg::hasSatisfiedTermConditions() const {
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


Real EMAlg::iterate( MaximizationStep &mstep ) {
    Real logZ = 0;
    Real likelihood = 0;

    _estep.run();
    logZ = _estep.logZ();

    // Expectation calculation
    for( Evidence::const_iterator e = _evidence.begin(); e != _evidence.end(); ++e ) {
        InfAlg* clamped = _estep.clone();
        e->applyEvidence( *clamped );
        clamped->init();
        clamped->run();

        likelihood += clamped->logZ() - logZ;

        mstep.addExpectations( *clamped );

        delete clamped;
    }

    // Maximization of parameters
    mstep.maximize( _estep.fg() );

    return likelihood;
}


Real EMAlg::iterate() {
    Real likelihood;
    for( size_t i = 0; i < _msteps.size(); ++i )
        likelihood = iterate( _msteps[i] );
    _lastLogZ.push_back( likelihood );
    ++_iters;
    return likelihood;
}


void EMAlg::run() {
    while( !hasSatisfiedTermConditions() )
        iterate();
}


} // end of namespace dai