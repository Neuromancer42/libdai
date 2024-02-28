#include <dai/dai_config.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <dai_ext/causal_bp.h>
#include <dai/util.h>
#include <dai/properties.h>

using namespace dai;

using namespace std;

//#define DAI_BP_FAST 1

static inline void CausalNormalize(Real &x, Real &y, bool logdomain = false) {
    if (logdomain) {
        Real e = y - x;
        DAI_ASSERT(!std::isnan(e));
        x = - dai::log(1 + dai::exp(e));
        y = - dai::log(1 + dai::exp(-e));
    } else {
        Real s = x + y;
        DAI_ASSERT(s > 0);
        x /= s;
        y /= s;
    }
}

static inline void CausalScale(Real &x, Real &y) {
    Real m = std::max(dai::abs(x), dai::abs(y));
    if (m != 0) {
        if (isinf(m)) {
            x = isinf(x) ? ((x < 0) ? -1 : 1) : 0;
            y = isinf(y) ? ((y < 0) ? -1 : 1) : 0;
        } else {
            x /= m;
            y /= m;
        }
    }
}

static inline void CausalScaleLog(Real &x, Real &y) {
    Real m = std::max(x, y);
    if (isinf(m)) {
        if (m > 0) {
            x = isinf(x) ? 0 : dai::log((Real) 0);
            y = isinf(y) ? 0 : dai::log((Real) 1);
        }
    } else {
        x -= m;
        y -= m;
    }
}
void CausalBP::setProperties(const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("tol") );
    DAI_ASSERT( opts.hasKey("logdomain") );
    DAI_ASSERT( opts.hasKey("updates") );

    props.tol = opts.getStringAs<Real>("tol");
    props.logdomain = opts.getStringAs<bool>("logdomain");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");

    if( opts.hasKey("maxiter") )
        props.maxiter = opts.getStringAs<size_t>("maxiter");
    else
        props.maxiter = 10000;
    if( opts.hasKey("maxtime") )
        props.maxtime = opts.getStringAs<Real>("maxtime");
    else
        props.maxtime = INFINITY;
    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<Real>("damping");
    else
        props.damping = 0.0;
    if( opts.hasKey("inference") )
        props.inference = opts.getStringAs<Properties::InfType>("inference");
    else
        props.inference = Properties::InfType::SUMPROD;
}

PropertySet CausalBP::getProperties() const {
    PropertySet opts;
    opts.set( "tol", props.tol );
    opts.set( "maxiter", props.maxiter );
    opts.set( "maxtime", props.maxtime );
    opts.set( "verbose", props.verbose );
    opts.set( "logdomain", props.logdomain );
    opts.set( "updates", props.updates );
    opts.set( "damping", props.damping );
    opts.set( "inference", props.inference );
    return opts;
}


string CausalBP::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "maxtime=" << props.maxtime << ",";
    s << "verbose=" << props.verbose << ",";
    s << "logdomain=" << props.logdomain << ",";
    s << "updates=" << props.updates << ",";
    s << "damping=" << props.damping << ",";
    s << "inference=" << props.inference << "]";
    return s.str();
}


void CausalBP::construct() {
    _edges.clear();
    _edges.reserve( nrVars() );
    _edge2lut.clear();
//    if( props.updates == Properties::UpdateType::SEQMAX )
//        _edge2lut.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i ) {
        _edges.push_back( { } );
        _edges[i].reserve( nbV(i).size() );
//        if( props.updates == Properties::UpdateType::SEQMAX ) {
//            _edge2lut.push_back( vector<LutType::iterator>() );
//            _edge2lut[i].reserve( nbV(i).size() );
//        }
        bforeach( const Neighbor &I, nbV(i) ) {
            EdgeProp newEP;
            newEP.message = Prob( var(i).states() );
            newEP.newMessage = Prob( var(i).states() );
            
//            if( DAI_BP_FAST ) {
//                newEP.index.reserve( factor(I).nrStates() );
//                for( IndexFor k( var(i), factor(I).vars() ); k.valid(); ++k )
//                    newEP.index.push_back( k );
//            }
            
            newEP.residual = 0.0;
            _edges[i].push_back( newEP );
//            if( props.updates == Properties::UpdateType::SEQMAX )
//                _edge2lut[i].push_back( _lut.insert( make_pair( newEP.residual, make_pair( i, _edges[i].size() - 1 ))) );
        }
    }

    _oldBeliefsV.clear();
    _oldBeliefsV.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i )
        _oldBeliefsV.push_back( Factor( var(i) ) );
//    _oldBeliefsF.clear();
//    _oldBeliefsF.reserve( nrFactors() );
//    for( size_t I = 0; I < nrFactors(); ++I )
//        _oldBeliefsF.push_back( Factor( factor(I).vars() ) );
    
    _updateSeq.clear();
    _updateSeq.reserve( nrEdges() );
    for( size_t I = 0; I < nrFactors(); I++ )
        bforeach( const Neighbor &i, nbF(I) )
            _updateSeq.push_back( Edge( i, i.dual ) );
}


void CausalBP::init() {
    for (size_t i = 0; i < nrVars(); ++i) {
        auto& vMsg = varMsgs[i];
        vMsg[0].reset(props.logdomain);
        vMsg[1].reset(props.logdomain);
    }
    Real c = props.logdomain ? 0.0 : 1.0;
    for( size_t i = 0; i < nrVars(); ++i ) {
        bforeach( const Neighbor &I, nbV(i) ) {
            message( i, I.iter ).fill( c );
            newMessage( i, I.iter ).fill( c );
//            if( props.updates == Properties::UpdateType::SEQMAX )
//                updateResidual( i, I.iter, 0.0 );
        }
    }
    _iters = 0;
}


void CausalBP::findMaxResidual( size_t &i, size_t &_I ) {
    DAI_ASSERT( !_lut.empty() );
    LutType::const_iterator largestEl = _lut.end();
    --largestEl;
    i  = largestEl->second.first;
    _I = largestEl->second.second;
}

//Prob CausalBP::calcIncomingMessageProduct(size_t I, bool without_i, size_t i ) const {
//    Factor Fprod( factor(I) );
//    Prob &prod = Fprod.p();
//    if( props.logdomain )
//        prod.takeLog();
//
//    // Calculate product of incoming messages and factor I
//    bforeach( const Neighbor &j, nbF(I) )
//        if( !(without_i && (j == i)) ) {
//            // prod_j will be the product of messages coming into j
//            Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 );
//            bforeach( const Neighbor &J, nbV(j) )
//                if( J != I ) { // for all J in nb(j) \ I
//                    if( props.logdomain )
//                        prod_j += message( j, J.iter );
//                    else
//                        prod_j *= message( j, J.iter );
//                }
//
//            // multiply prod with prod_j
////            if( !DAI_BP_FAST ) {
//                // UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION
//                if( props.logdomain )
//                    Fprod += Factor( var(j), prod_j );
//                else
//                    Fprod *= Factor( var(j), prod_j );
////            } else {
////                // OPTIMIZED VERSION
////                size_t _I = j.dual;
////                // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
////                const ind_t &ind = index(j, _I);
////                
////                for( size_t r = 0; r < prod.size(); ++r )
////                    if( props.logdomain )
////                        prod.set( r, prod[r] + prod_j[ind[r]] );
////                    else
////                        prod.set( r, prod[r] * prod_j[ind[r]] );
////            }
//        }
//    return prod;
//}


void CausalBP::calcNewMessage(size_t i, size_t _I ) {
    size_t I = nbV(i, _I);
//    Prob marg;
    Real marg0, marg1;
    switch (factor(I).type) {
        case CausalFactor::Singleton: {
//            Real p = factor(I).prob();
//            marg = Prob(std::vector<Real>{1.0 - p, p});
            marg1 = factor(I).prob();
            marg0 = static_cast<Real>(1.0) - marg1;
            break;
        }
        case CausalFactor::DefiniteAnd: {
//            Prob mask(var(i).states(), 1.0 );
            Real mask0 = 1.0, mask1 = 1.0;
            if (factor(I).head_clamped) {
//                mask *= factor(I).head_mask.exp(1 / static_cast<Real>(factor(I).body().size()));
                mask0 = factor(I).head_mask[0];
                mask1 = factor(I).head_mask[1];
            }
            if (factor(I).head() == var(i)) {
                Real t0 = 1, t1 = 1, e1 = 0;

                bforeach(const Neighbor &j, nbF(I)) {
                    if (j != i) {
//                        Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0);
//                        DAI_ASSERT(prod_j.size() == 2);
                        auto &vMsg = varMsgs[j];
                        auto &origMsg = message(j, j.dual);
                        Real prod_j0 = vMsg[0].residual(props.logdomain, I, origMsg[0]);
                        Real prod_j1 = vMsg[1].residual(props.logdomain, I, origMsg[1]);
                        
                        if (props.logdomain) {
//                            prod_j.takeExp();
                            prod_j0 = dai::exp(prod_j0);
                            prod_j1 = dai::exp(prod_j1);
//                            prod_j.normalize();
                        }
                        CausalNormalize(prod_j0, prod_j1);

                        constexpr Real arg0 = 0, arg1 = 1;
//                        Real a0 = (prod_j[0] + prod_j[1]);
                        Real a0 = prod_j0 + prod_j1;
//                        Real a1 = (arg0 * prod_j[0] + arg1 * prod_j[1]);
                        Real a1 = arg0 * prod_j0 + arg1 * prod_j1;
//                        Real delta = (1-arg0)*prod_j[0] + (1-arg1)*prod_j[1];
                        Real delta = (1 - arg0) * prod_j0 + (1  - arg1) * prod_j1;
                        t0 *= a0;
                        t1 *= a1;
                        CausalScale(t0, t1);
                        if (a1 != 0 && a0 == a1 && delta != 0) {
                            e1 += delta / a1;
                        }
                    }
                }
//                marg = Prob(std::vector<Real>{e1 * t0 + (t0 - t1), t1}) * mask;
                marg0 = (e1 * t0 + (t0 - t1)) * mask0;
                marg1 = t1 * mask1;
            } else {
                Real t0 = 1, t1 = 1;
                bforeach(const Neighbor &j, nbF(I)) {
                    if (j != i) {
//                        Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0);
//                        DAI_ASSERT(prod_j.size() == 2);
                        auto &vMsg = varMsgs[j];
                        auto &origMsg = message(j, j.dual);
                        Real prod_j0 = vMsg[0].residual(props.logdomain, I, origMsg[0]);
                        Real prod_j1 = vMsg[1].residual(props.logdomain, I, origMsg[1]);
                        
                        if (props.logdomain) {
//                            prod_j.takeExp();
                            prod_j0 = dai::exp(prod_j0);
                            prod_j1 = dai::exp(prod_j1);
//                            prod_j.normalize();
                        }
                        CausalNormalize(prod_j0, prod_j1);
                        
                        if (factor(I).head() == var(j)) {
                            t0 *= (prod_j1 * mask1 - prod_j0 * mask0);
                            t1 *= prod_j0 * mask0;
                        } else {
                            constexpr Real arg0 = 0, arg1 = 1;
                            t0 *= (arg0 * prod_j0 + arg1 * prod_j1);
                            t1 *= prod_j0 + prod_j1;
                        }
                        CausalScale(t0, t1);
                    }
                }
                Real px0, px1;
                constexpr Real arg0 = 0, arg1 = 1;
                px0 = arg0 * t0 + t1;
                px1 = arg1 * t0 + t1;
//                marg = Prob(std::vector<Real>{px0, px1});
                marg0 = px0;
                marg1 = px1;
            }
            break;
        }
        case CausalFactor::DefiniteOr: {
//            Prob mask(var(i).states(), 1.0 );
            Real mask0 = 1.0, mask1 = 1.0;
            if (factor(I).head_clamped) {
//                mask *= factor(I).head_mask;
                mask0 = factor(I).head_mask[0];
                mask1 = factor(I).head_mask[1];
            }
            if (factor(I).head() == var(i)) {
                Real t0 = 1, t1 = 1, e1 = 0;
                bforeach(const Neighbor &j, nbF(I)) {
                    if (j != i) {
//                        Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 );
//                        DAI_ASSERT(prod_j.size() == 2);
                        auto &vMsg = varMsgs[j];
                        auto &origMsg = message(j, j.dual);
                        Real prod_j0 = vMsg[0].residual(props.logdomain, I, origMsg[0]);
                        Real prod_j1 = vMsg[1].residual(props.logdomain, I, origMsg[1]);

                        if (props.logdomain) {
//                            prod_j.takeExp();
                            prod_j0 = dai::exp(prod_j0);
                            prod_j1 = dai::exp(prod_j1);
//                            prod_j.normalize();
                        }
                        CausalNormalize(prod_j0, prod_j1);
                        
                        constexpr Real arg0 = 1, arg1 = 0;
                        Real a0 = prod_j0 + prod_j1;
                        Real a1 = arg0 * prod_j0 + arg1 * prod_j1;
                        Real delta = (1 - arg0) * prod_j0 + (1 - arg1) * prod_j1;
                        t0 *= a0;
                        t1 *= a1;
                        CausalScale(t0, t1);
                        if (a1 != 0 && a0 == a1 && delta != 0) {
                            e1 += delta / a1;
                        }
                    }
                }

//                marg = Prob(std::vector<Real>{t1, e1 * t0 + (t0 - t1)}) * mask;
                marg0 = t1 * mask0;
                marg1 = (e1 * t0 + (t0 - t1)) * mask1;
            } else {
                Real t0 = 1, t1 = 1;
                bforeach(const Neighbor &j, nbF(I)) {
                    if (j != i) {
//                        Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 );
//                        DAI_ASSERT(prod_j.size() == 2);
                        auto &vMsg = varMsgs[j];
                        auto &origMsg = message(j, j.dual);
                        Real prod_j0 = vMsg[0].residual(props.logdomain, I, origMsg[0]);
                        Real prod_j1 = vMsg[1].residual(props.logdomain, I, origMsg[1]);

                        if (props.logdomain) {
//                            prod_j.takeExp();
                            prod_j0 = dai::exp(prod_j0);
                            prod_j1 = dai::exp(prod_j1);
//                            prod_j.normalize();
                        }
                        CausalNormalize(prod_j0, prod_j1);
                        
                        if (factor(I).head() == var(j)) {
                            t0 *= prod_j0 * mask0 - prod_j1 * mask1;
                            t1 *= prod_j1 * mask1;
                        } else {
                            constexpr Real arg0 = 1, arg1 = 0;
                            t0 *= arg0 * prod_j0 + arg1 * prod_j1;
                            t1 *= prod_j0 + prod_j1;
                        }
                        CausalScale(t0, t1);
                    }
                }
                Real px0, px1;
                constexpr Real arg0 = 1, arg1 = 0;
                px0 = arg0 * t0 + t1;
                px1 = arg1 * t0 + t1;
//                marg = Prob(std::vector<Real>{px0, px1});
                marg0 = px0;
                marg1 = px1;
            }
            break;
        }
    }
    
//    if( factor(I).vars().size() == 1 ) 
//        marg = factor(I).p();
//    else {
//        Factor Fprod(factor(I));
//        Prob &prod = Fprod.p();
//        prod = calcIncomingMessageProduct(I, true, i);
//
//        if (props.logdomain) {
//            prod -= prod.max();
//            prod.takeExp();
//        }
//
//        if (!DAI_BP_FAST) {
//            if (props.inference == Properties::InfType::SUMPROD)
//                marg = Fprod.marginal(var(i)).p();
//            else
//                marg = Fprod.maxMarginal(var(i)).p();
//        } else {
//            marg = Prob(var(i).states(), 0.0);
//            const ind_t ind = index(i, _I);
//            if (props.inference == Properties::InfType::SUMPROD)
//                for (size_t r = 0; r < prod.size(); ++r)
//                    marg.set(ind[r], marg[ind[r]] + prod[r]);
//            else
//                for (size_t r = 0; r < prod.size(); ++r)
//                    if (prod[r] > marg[ind[r]])
//                        marg.set(ind[r], prod[r]);
//            marg.normalize();
//        }
//    }
        
//    marg.normalize();
    CausalNormalize(marg0, marg1);
    if( props.logdomain ) {
        newMessage(i,_I).set(0, dai::log(marg0));
        newMessage(i,_I).set(1, dai::log(marg1));
    } else {
        newMessage(i,_I).set(0, marg0);
        newMessage(i,_I).set(1, marg1);
    }
//    std::cerr << "newMessage(" << i << "<-" << I << ")=" << marg << std::endl;
    
//    if( props.updates == Properties::UpdateType::SEQMAX )
//        updateResidual( i, _I , dist( newMessage( i, _I ), message( i, _I ), DISTLINF ) );
}


Real CausalBP::run() {
    if( props.verbose >= 1 )
        cerr << "Starting " << identify() << "...";
    if( props.verbose >= 3)
        cerr << endl;
    
    double tic = toc();
    
    Real maxDiff = INFINITY;
    for( ; _iters < props.maxiter && maxDiff > props.tol && (toc() - tic) < props.maxtime; _iters++ ) {
        /* if( props.updates == Properties::UpdateType::SEQMAX ) {
            if( _iters == 0 ) {
                for( size_t i = 0; i < nrVars(); ++i )
                    bforeach( const Neighbor &I, nbV(i) )
                        calcNewMessage( i, I.iter );
            }
            for( size_t t = 0; t < _updateSeq.size(); ++t ) {
                size_t i, _I;
                findMaxResidual( i, _I );
                updateMessage( i, _I );

                bforeach( const Neighbor &J, nbV(i) ) {
                    if( J.iter != _I ) {
                        bforeach( const Neighbor &j, nbF(J) ) {
                            size_t _J = j.dual;
                            if( j != i )
                                calcNewMessage( j, _J );
                        }
                    }
                }
            }
        } else */ if( props.updates == Properties::UpdateType::PARALL ) {
            for( size_t i = 0; i < nrVars(); ++i )
                bforeach( const Neighbor &I, nbV(i) )
                    calcNewMessage( i, I.iter );

            for( size_t i = 0; i < nrVars(); ++i )
                bforeach( const Neighbor &I, nbV(i) )
                    updateMessage( i, I.iter );

            for (size_t i = 0; i < nrVars(); ++i) {
                auto & vMsg = varMsgs[i];
                vMsg[0].reset(props.logdomain);
                vMsg[1].reset(props.logdomain);
                bforeach( const Neighbor &I, nbV(i) ) {
                    auto & m = message(i, I.iter);
                    vMsg[0].accumulate(props.logdomain, I, m[0]);
                    vMsg[1].accumulate(props.logdomain, I, m[1]);
                    if (props.logdomain) {
                        CausalScaleLog(vMsg[0].msg, vMsg[1].msg);
                    } else {
                        CausalScale(vMsg[0].msg, vMsg[1].msg);
                    }
                }
            }
        } /* else {
            if( props.updates == Properties::UpdateType::SEQRND )
                // TOOD: modify randome_shuffle to support multi-thread call
                random_shuffle( _updateSeq.begin(), _updateSeq.end(), rnd );
            
            bforeach( const Edge &e, _updateSeq ) {
                calcNewMessage( e.first, e.second );
                updateMessage( e.first, e.second );
            }
        } */
        
        maxDiff = -INFINITY;
        for( size_t i = 0; i < nrVars(); ++i ) {
            Factor b( beliefV(i) );
            maxDiff = std::max( maxDiff, dist( b, _oldBeliefsV[i], DISTLINF ) );
            _oldBeliefsV[i] = b;
        }
//        for( size_t I = 0; I < nrFactors(); ++I ) {
//            Factor b( beliefF(I) );
//            // TODO: how to calculate this?
//            maxDiff = std::max( maxDiff, dist( b, _oldBeliefsF[I], DISTLINF ) );
//            _oldBeliefsF[I] = b;
//        }
        
        if( props.verbose >= 3 )
            cerr << name() << "::run:  maxdiff " << maxDiff << " after " << _iters+1 << " passes" << endl;
    }
    
    if( maxDiff > _maxdiff )
        _maxdiff = maxDiff;
    
    if( props.verbose >= 1 ) {
        if( maxDiff > props.tol ) {
            if( props.verbose == 1 )
                cerr << endl;
                cerr << name() << "::run:  WARNING: not converged after " << _iters << " passes (" << toc() - tic << " seconds)...final maxdiff:" << maxDiff << endl;
        } else {
            if( props.verbose >= 3 )
                cerr << name() << "::run:  ";
                cerr << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }
    
    return maxDiff;
}


void CausalBP::calcBeliefV(size_t i, Prob &p ) const {
    p = Prob( var(i).states(), props.logdomain ? 0.0 : 1.0 );
    bforeach( const Neighbor &I, nbV(i) )
        if( props.logdomain )
            p += newMessage( i, I.iter );
        else
            p *= newMessage( i, I.iter );
}

Factor CausalBP::beliefV(size_t i ) const {
    Prob p;
    calcBeliefV( i, p );
    
    if( props.logdomain ) {
        p -= p.max();
        p.takeExp();
    }
    p.normalize();
    
    return( Factor( var(i), p ) );
}

//Factor CausalBP::beliefF(size_t I ) const {
//    Prob p;
//    calcBeliefF( I, p );
//
//    if( props.logdomain ) {
//        p -= p.max();
//        p.takeExp();
//    }
//    p.normalize();
//
//    return( Factor( factor(I).vars(), p ) );
//}


vector<Factor> CausalBP::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); ++i )
        result.push_back( beliefV(i) );
//    for( size_t I = 0; I < nrFactors(); ++I )
//        result.push_back( beliefF(I) );
    return result;
}


Factor CausalBP::belief(const VarSet &ns ) const {
    if( ns.size() == 0 )
        return Factor();
    else if( ns.size() == 1 )
        return beliefV( findVar( *(ns.begin() ) ) );
    else {
        DAI_THROW(BELIEF_NOT_AVAILABLE);

//        size_t I;
//        for( I = 0; I < nrFactors(); I++ )
//            if( factor(I).vars() >> ns )
//                break;
//        if( I == nrFactors() )
//            DAI_THROW(BELIEF_NOT_AVAILABLE);
//        return beliefF(I).marginal(ns);
    }
}


Real CausalBP::logZ() const {
    Real sum = 0.0;
    for( size_t i = 0; i < nrVars(); ++i )
        sum += (1.0 - nbV(i).size()) * beliefV(i).entropy();
//    for( size_t I = 0; I < nrFactors(); ++I )
//        sum -= calcFactorDistKL(I);
    return sum;
}


void CausalBP::init(const VarSet &ns ) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); ++n ) {
        size_t ni = findVar( *n );
        auto & vMsg = varMsgs[ni];
        vMsg[0].reset(props.logdomain);
        vMsg[1].reset(props.logdomain);
        bforeach( const Neighbor &I, nbV( ni ) ) {
            Real val = props.logdomain ? 0.0 : 1.0;
            message( ni, I.iter ).fill( val );
            newMessage( ni, I.iter ).fill( val );
//            if( props.updates == Properties::UpdateType::SEQMAX )
//                updateResidual( ni, I.iter, 0.0 );
        }
    }
    _iters = 0;
}


void CausalBP::updateMessage(size_t i, size_t _I ) {
    if( recordSentMessages )
        _sentMessages.push_back(make_pair(i,_I));
    if( props.damping == 0.0 ) {
        message(i,_I) = newMessage(i,_I);
//        if( props.updates == Properties::UpdateType::SEQMAX )
//            updateResidual( i, _I, 0.0 );
    } else {
        if( props.logdomain )
            message(i,_I) = (message(i,_I) * props.damping) + (newMessage(i,_I) * (1.0 - props.damping));
        else
            message(i,_I) = (message(i,_I) ^ props.damping) * (newMessage(i,_I) ^ (1.0 - props.damping));
//        if( props.updates == Properties::UpdateType::SEQMAX )
//            updateResidual( i, _I, dist( newMessage(i,_I), message(i,_I), DISTLINF ) );
    }
}


void CausalBP::updateResidual(size_t i, size_t _I, Real r ) {
    EdgeProp* pEdge = &_edges[i][_I];
    pEdge->residual = r;

    // rearrange look-up table (delete and reinsert new key)
    _lut.erase( _edge2lut[i][_I] );
    _edge2lut[i][_I] = _lut.insert( make_pair( r, make_pair(i, _I) ) );
}

//Real CausalBP::calcFactorDistKL(size_t I) {
//    Real dist;
//    switch (factor(I).type) {
//        case CausalFactor::Singleton: {
//            break;
//        }
//        case CausalFactor::DefiniteAnd: {
//            break;
//        }
//        case CausalFactor::DefiniteOr: {
//            break;
//        }
//    } 
//    return dist;
//}
