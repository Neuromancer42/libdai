//
// Created by Yifan Chen on 2024/1/1.
//

#ifndef LIBDAI_FACTORIZED_FG_H
#define LIBDAI_FACTORIZED_FG_H


#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <dai/prob.h>
#include <dai/varset.h>
#include <dai/index.h>
#include <dai/util.h>
#include <dai/bipgraph.h>
#include <dai/graph.h>

namespace dai {
    
    class CausalFactor {
    public:
        enum CausalType {
            Singleton = 'I',
            DefiniteAnd = '*',
            DefiniteOr = '+'
        };
        
        CausalType type;

        bool head_clamped;
        Prob head_mask;
        
    private:
        Var _head;
        VarSet _body; 
        Real _p;
        VarSet _vs; // duplicate for usage
        
        friend std::ostream& operator<< (std::ostream&, const CausalFactor&);
        
    public:
        explicit CausalFactor( const Var &v, Real p = 1 ) : _head(v), _body(), _p(p), type(Singleton) {
            _vs = VarSet();
            _vs.insert(_head);
            head_clamped = false;
            head_mask = Prob{std::vector<Real>{1, 1}};
        }
        
        CausalFactor( const Var &head, const VarSet &body, bool isAnd, Real p = 1) : _head(head), _body(body), _p(p) {
            if (isAnd)
                type = DefiniteAnd;
            else
                type = DefiniteOr;
            _vs = VarSet(_body);
            _vs.insert(_head);
            head_clamped = false;
        }
        
        CausalFactor gen_clamped(Var i, size_t x) const {
            CausalFactor newFac(*this);
            switch (type) {
                case Singleton:
                    if (i == _head)
                        newFac._p = static_cast<Real>(x);
                    else {
                        std::cerr << "Clamp unrelated var " << i << " in factor " << *this << std::endl;
                    }
                    break;
                case DefiniteAnd:
                case DefiniteOr:
                    if (i == _head) {
                        newFac.head_clamped = true;
                        newFac.head_mask = Prob(2, static_cast<Real>(0));
                        newFac.head_mask.set(x, static_cast<Real>(1));
                        if (head_clamped) {
                            std::cerr << "Duplicate clamp head of factor " << *this << "with val " << x << std::endl;
                            newFac.head_mask *= head_mask;
                        }
                    } else {
                        if (_body.contains(i)) {
                            newFac._body.erase(i);
                        }
                    }
                    break;
            }
            return newFac;
        }
        
        Real prob() const { return _p; }
        
        const Var& head() const { return _head; }
        const VarSet& body() const { return _body; }
        
        const VarSet& vars() const {
            return _vs;
        }

        std::string toString() const {
            std::ostringstream ss;
            ss << *this;
            return ss.str();
        }
    };
    
    class CausalFactorGraph {
    private:
        BipartiteGraph _G;
        std::vector<Var> _vars;
        std::vector<CausalFactor> _factors;
        std::map<size_t, CausalFactor> _backup;
        
    public:
        CausalFactorGraph() : _G(), _vars(), _factors(), _backup() {}
        
        CausalFactorGraph( const std::vector<CausalFactor>& P );

        template<typename CausalFactorInputIterator, typename VarInputIterator>
        CausalFactorGraph(CausalFactorInputIterator facBegin, CausalFactorInputIterator facEnd, VarInputIterator varBegin, VarInputIterator varEnd, size_t nrFacHint = 0, size_t nrVarHint = 0 );
        
        virtual ~CausalFactorGraph() {}
        
        virtual CausalFactorGraph* clone() const { return new CausalFactorGraph(*this); }
        
        const Var& var( size_t i ) const {
            DAI_DEBASSERT( i < nrVarrs() );
            return _vars[i];
        }
        const std::vector<Var>& vars() const { return _vars; }
        
        const CausalFactor& factor( size_t I ) const {
            DAI_DEBASSERT( I < nrFactors() );
            return _factors[I];
        }
        const std::vector<CausalFactor>& factors() const { return _factors; }
        const Neighbors& nbV( size_t i ) const { return _G.nb1(i); }
        const Neighbors& nbF( size_t I ) const { return _G.nb2(I); }
        const Neighbor& nbV( size_t i, size_t _I ) const { return _G.nb1(i)[_I]; }
        const Neighbor& nbF( size_t I, size_t _i ) const { return _G.nb2(I)[_i]; }
        
        const BipartiteGraph& bipGraph() const { return _G; }
        
        size_t nrVars() const { return vars().size(); }
        size_t nrFactors() const { return factors().size(); }
        size_t nrEdges() const { return _G.nrEdges(); }
        
        size_t findVar( const Var& n ) const {
            size_t i = find( vars().begin(), vars().end(), n ) - vars().begin();
            if( i == nrVars() )
                DAI_THROW(OBJECT_NOT_FOUND);
            return i;
        }
        
        SmallSet<size_t> findVars( const VarSet& ns ) const {
            SmallSet<size_t> result;
            for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
                result.insert( findVar( *n ) );
            return result;
        }
        
        size_t findFactor( const VarSet& ns ) const {
            size_t I;
            for( I = 0; I < nrFactors(); I++ )
                if( factor(I).vars() == ns )
                    break;
            if( I == nrFactors() )
                DAI_THROW(OBJECT_NOT_FOUND);
            return I;
        }
        
        VarSet inds2vars( const std::vector<size_t>& inds ) {
            VarSet vs;
            for( std::vector<size_t>::const_iterator it = inds.begin(); it != inds.end(); it++ )
                vs.insert( var(*it) );
            return vs;
        }

        VarSet Delta( size_t i ) const;
        VarSet Delta( const VarSet& vs ) const;
        
        VarSet delta( size_t i ) const {
            return( Delta( i ) / var( i ) );
        }
        VarSet delta( const VarSet& vs ) const {
            return Delta( vs ) / vs;
        }
        
        SmallSet<size_t> Deltai( size_t i ) const;
        SmallSet<size_t> deltai( size_t i ) const {
            return( Deltai( i ) / i );
        }
        
        bool isConnected() const { return _G.isConnected(); }
        bool isTree() const { return _G.isTree(); }
        bool isPairwise() const;
        bool isBinary() const { return true; };

        GraphAL MarkovGraph() const;
        
        bool isMaximal( size_t I ) const;
        size_t maximalFactor( size_t I ) const;
        std::vector<VarSet> maximalFactorDomains() const;
        Real logScore( const std::vector<size_t>& statevec ) const;

        virtual void setFactor( size_t I, const CausalFactor& newFactor, bool backup = false ) {
            DAI_ASSERT( newFactor.vars() << factor(I).vars() );
            if( backup )
                backupFactor( I );
            _factors[I] = newFactor;
        }

        virtual void setFactors( const std::map<size_t, CausalFactor>& facs, bool backup = false ) {
            for( std::map<size_t, CausalFactor>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ ) {
                if( backup )
                    backupFactor( fac->first );
                setFactor( fac->first, fac->second );
            }
        }

        void backupFactor( size_t I );
        void restoreFactor( size_t I );

        virtual void backupFactors( const std::set<size_t>& facs );
        virtual void restoreFactors();

        void backupFactors( const VarSet& ns );
        void restoreFactors( const VarSet& ns );

        CausalFactorGraph maximalFactors() const;

        CausalFactorGraph clamped( size_t i, size_t x ) const;
        
        virtual void clamp( size_t i, size_t x, bool backup = false );
        
        void clampVar( size_t i, const std::vector<size_t>& xis, bool backup = false );
        
        void clampFactor( size_t I, const std::vector<size_t>& xIs, bool backup = false );

        virtual void makeCavity( size_t i, bool backup = false );
        virtual void makeRegionCavity( std::vector<size_t> facInds, bool backup );
        
        virtual void ReadFromFile( const char *filename );
        virtual void WriteToFile( const char *filename, size_t precision=15 ) const;
        
        friend std::ostream& operator<< (std::ostream& os, const CausalFactorGraph& fg );
        friend std::istream& operator>> (std::istream& is, CausalFactorGraph& fg );

        virtual void printDot( std::ostream& os ) const;
        
        std::string toString() const {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }

        void fromString( const std::string& s ) {
            std::stringstream ss( s );
            ss >> *this;
        }

    private:
        /// Part of constructors (creates edges, neighbors and adjacency matrix)
        void constructGraph( size_t nrEdges );
    };

    template<typename CausalFactorInputIterator, typename VarInputIterator>
    CausalFactorGraph::CausalFactorGraph(CausalFactorInputIterator facBegin, CausalFactorInputIterator facEnd, VarInputIterator varBegin, VarInputIterator varEnd, size_t nrFacHint, size_t nrVarHint ) : _G(), _backup() {
        // add factors
        size_t nrEdges = 0;
        _factors.reserve( nrFacHint );
        for( CausalFactorInputIterator p2 = facBegin; p2 != facEnd; ++p2 ) {
            _factors.push_back( *p2 );
            nrEdges += p2->vars().size();
        }

        // add variables
        _vars.reserve( nrVarHint );
        for( VarInputIterator p1 = varBegin; p1 != varEnd; ++p1 )
            _vars.push_back( *p1 );

        // create graph structure
        constructGraph( nrEdges );
    }
}
#endif //LIBDAI_FACTORIZED_FG_H
