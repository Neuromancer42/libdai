//
// Created by Yifan Chen on 2023/12/24.
//

#ifndef LIBDAI_CAUSAL_BP_H
#define LIBDAI_CAUSAL_BP_H

#include "dai/alldai.h"
#include "dai_ext/causal_fg.h"

namespace dai {
    class CausalBP : public DAIAlg<CausalFactorGraph> {
    protected:
//        typedef std::vector<size_t> ind_t;
        struct EdgeProp {
            Prob message;
            Prob newMessage;
            Real residual;
        };
        
        std::vector<std::vector<EdgeProp> > _edges;
        
//        typedef std::multimap<Real, std::pair<size_t, size_t> > LutType;
//        std::vector<std::vector<LutType::iterator> > _edge2lut;
//        LutType _lut;

        Real _maxdiff;
        size_t _iters;
//        std::vector<std::pair<size_t, size_t> > _sentMessages;
        std::vector<Factor> _oldBeliefsV;
//        std::vector<Factor> _oldBeliefsF;
//        std::vector<Edge> _updateSeq;

        std::vector<double> _lowPassBeliefs;

    public:
        /// Parameters are the same as BP
        struct Properties {
            /// Enumeration of possible update schedules
            /** The following update schedules have been defined:
             *  - PARALL parallel updates
             *  - SEQFIX sequential updates using a fixed sequence
             *  - SEQRND sequential updates using a random sequence
             *  - SEQMAX maximum-residual updates [\ref EMK06]
             */
            DAI_ENUM(UpdateType/*,SEQFIX,SEQRND,SEQMAX*/,PARALL);

            /// Enumeration of inference variants
            /** There are two inference variants:
             *  - SUMPROD Sum-Product
             *  - MAXPROD Max-Product (equivalent to Min-Sum)
             */
            DAI_ENUM(InfType,SUMPROD/*,MAXPROD*/);

            /// Verbosity (amount of output sent to stderr)
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Maximum time (in seconds)
            double maxtime;

            /// Tolerance for convergence test
            Real tol;

            /// Whether updates should be done in logarithmic domain or not
            bool logdomain;

            /// Damping constant (0.0 means no damping, 1.0 is maximum damping)
            Real damping;

            /// Message update schedule
            UpdateType updates;

            /// Inference variant
            InfType inference;
        } props;

//        /// Specifies whether the history of message updates should be recorded
//        bool recordSentMessages;
        
    public:
        CausalBP() : DAIAlg<CausalFactorGraph>(),  
                _edges()/*, _edge2lut(), _lut()*/, _maxdiff(0.0), _iters(0U)/*, recordSentMessages(false),
                _sentMessages()*/, _oldBeliefsV()/*, _oldBeliefsF()*//*, _updateSeq()*/, props(), 
//                factorMsgs(nrFactors(), {AccumulateMsg{props.logdomain}, AccumulateMsg{props.logdomain}}),
                varMsgs(nrVars(), {AccumulateMsg{props.logdomain}, AccumulateMsg{props.logdomain}}) {}
        
        CausalBP(const CausalFactorGraph &x, const PropertySet &opts ) : DAIAlg<CausalFactorGraph>(x),
                _edges()/*, _edge2lut(), _lut()*/, _maxdiff(0.0), _iters(0U)/*, recordSentMessages(false),
                _sentMessages()*/, _oldBeliefsV()/*, _oldBeliefsF()*//*, _updateSeq()*/, props(), 
//                factorMsgs(nrFactors(), {AccumulateMsg{props.logdomain}, AccumulateMsg{props.logdomain}}),
                varMsgs(nrVars(), {AccumulateMsg{props.logdomain}, AccumulateMsg{props.logdomain}}) {
            setProperties( opts );
            construct();
        }
        
        CausalBP( const CausalBP &x ) : DAIAlg<CausalFactorGraph>(x), 
                _edges(x._edges)/*, _edge2lut(x._edge2lut), _lut(x._lut)*/, _maxdiff(x._maxdiff), _iters(x._iters)/*, recordSentMessages(x.recordSentMessages), 
                _sentMessages(x._sentMessages)*/, _oldBeliefsV(x._oldBeliefsV)/*, _updateSeq(x._updateSeq)*/, props(x.props), 
                _lowPassBeliefs(x._lowPassBeliefs), varMsgs(x.varMsgs) {
//            for( LutType::iterator l = _lut.begin(); l != _lut.end(); ++l )
//                _edge2lut[l->second.first][l->second.second] = l;
        }
        
        CausalBP& operator=(const CausalBP &x ) {
            if (this != &x) {
                InfAlg::operator=( x );
                _edges = x._edges;
//                _lut = x._lut;
//                for( LutType::iterator l = _lut.begin(); l != _lut.end(); ++l )
//                    _edge2lut[l->second.first][l->second.second] = l;
                _maxdiff = x._maxdiff;
                _iters = x._iters;
//                recordSentMessages = x.recordSentMessages;
//                _sentMessages = x._sentMessages;
                _oldBeliefsV = x._oldBeliefsV;
//                _updateSeq = x._updateSeq;
                props = x.props;
                _lowPassBeliefs = x._lowPassBeliefs;
                varMsgs = x.varMsgs;
            }
            return *this;
        }
        
        CausalBP* clone() const override { return new CausalBP(*this); }
        CausalBP* construct(const FactorGraph &fg, const PropertySet &opts ) const override {
            DAI_THROW(NOT_IMPLEMENTED);
            return nullptr;
        }
       
        std::string name() const override { return "CausalBP"; }
        Factor belief( const Var &v ) const override { return beliefV( findVar(v) ); }
        Factor belief( const VarSet &vs ) const override;
        Factor beliefV( size_t i ) const override;
        Factor beliefF( size_t I ) const override;
        std::vector<Factor> beliefs() const override;
        Real logZ() const override;
        
        std::vector<size_t> findMaximum() const override { return dai::findMaximum( *this ); }
        void init() override;
        void init( const VarSet &ns ) override;
        Real run() override;
        Real maxDiff() const override { return _maxdiff; }
        size_t Iterations() const override { return _iters; }
        void setMaxIter( size_t maxiter ) override { props.maxiter = maxiter; }
        void setProperties( const PropertySet &opts ) override;
        PropertySet getProperties() const override;
        std::string printProperties() const override;
        
//        const std::vector<std::pair<size_t, size_t> >& getSentMessages() const {
//            return _sentMessages;
//        }
//        
//        void clearSentMessages() { _sentMessages.clear(); }
        

        double run(double tolerance, size_t minIters, size_t maxIters, size_t histLength);
        double newBelief(size_t varIndex) const {
            return _lowPassBeliefs[varIndex];
        }
    protected:
        const Prob & message(size_t i, size_t _I) const { return _edges[i][_I].message; }
        Prob & message(size_t i, size_t _I) { return _edges[i][_I].message; }
        const Prob & newMessage(size_t i, size_t _I) const { return _edges[i][_I].newMessage; }
        Prob & newMessage(size_t i, size_t _I) { return _edges[i][_I].newMessage; }
        const Real & residual(size_t i, size_t _I) const { return _edges[i][_I].residual; }
        Real & residual(size_t i, size_t _I) { return _edges[i][_I].residual; }
        
//        virtual Prob calcIncomingMessageProduct( size_t I, bool withoud_i, size_t i ) const;
        virtual void calcNewMessage( size_t i, size_t _I );
        void updateMessage( size_t i, size_t _I );
        void updateResidual( size_t i, size_t _I, Real r );
        void findMaxResidual( size_t &i, size_t &_I );
        virtual void calcBeliefV( size_t i, Prob &p ) const;
//        virtual void calcBeliefF( size_t I, Prob &p ) const {
//            p = calcIncomingMessageProduct( I, false, 0 );
//        }
        
        virtual void construct();
        
    private:
//        Real calcFactorDistKL(size_t I);
        
        struct AccumulateMsg {
            Real msg;
            std::set<size_t> zeros{};
            explicit AccumulateMsg(bool logdomain) {
                msg = logdomain ? 0.0 : 1.0;
            }
            void reset(size_t id, Real origMsg, bool logdomain) {
                if (zeros.erase(id) == 0) {
                    if (logdomain)
                        msg -= origMsg;
                    else
                        msg /= origMsg;
                }
            }
            void reset(bool logdomain) {
                msg = logdomain ? 0.0 : 1.0;
                zeros.clear();
            }
            void accumulate(bool logdomain, size_t id, Real m) {
                if (logdomain) {
                    if (isinf(m))
                        zeros.insert(id);
                    else
                        msg += m;
                } else {
                    if (m == 0)
                        zeros.insert(id);
                    else
                        msg *= m;
                }
            }
            Real residual(bool logdomain, size_t i, Real m) const {
                if (zeros.size() > 1)
                    return logdomain ? dai::log((Real) 0) : 0;
                if (zeros.size() == 1) {
                    if (zeros.find(i) != zeros.end())
                        return msg;
                    return logdomain ? dai::log((Real) 0) : 0;
                }
                return logdomain ? (msg - m) : (msg / m);
            }
            Real get(bool logdomain) const {
                if (zeros.empty())
                    return msg;
                return logdomain ? dai::log((Real) 0) : 0;
            }
        };
        std::vector<std::array<AccumulateMsg, 2>> /*factorMsg,*/ varMsgs;
    };
}
#endif //LIBDAI_CAUSAL_BP_H
