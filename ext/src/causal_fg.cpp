//
// Created by Yifan Chen on 2024/1/10.
//
#include <iostream>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>
#include <dai_ext/causal_fg.h>
#include <dai/util.h>
#include <dai/exceptions.h>
#include <boost/lexical_cast.hpp>

using namespace dai;

using namespace std;

std::ostream& dai::operator<< (std::ostream& os, const CausalFactor& f) {
    os << "(" << f.head() << " " << static_cast<char>(f.type);
    if (f.type == CausalFactor::Singleton) {
        os << " " << std::setw(os.precision()+4) << f.prob();
    } else {
        for (const auto & v : f.body()) {
            os << " " << v.label();
        }
    }
    os << ")";
    return os;
}

CausalFactorGraph::CausalFactorGraph( const std::vector<CausalFactor> &P ) : _G(), _backup() {
    // add factors, obtain variables
    set<Var> varset;
    _factors.reserve( P.size() );
    size_t nrEdges = 0;
    for( vector<CausalFactor>::const_iterator p2 = P.begin(); p2 != P.end(); p2++ ) {
        _factors.push_back( *p2 );
        copy( p2->vars().begin(), p2->vars().end(), inserter( varset, varset.begin() ) );
        nrEdges += p2->vars().size();
    }

    // add vars
    _vars.reserve( varset.size() );
    for( set<Var>::const_iterator p1 = varset.begin(); p1 != varset.end(); p1++ )
        _vars.push_back( *p1 );

    // create graph structure
    constructGraph( nrEdges );
}


void CausalFactorGraph::constructGraph( size_t nrEdges ) {
    // create a mapping for indices
    hash_map<size_t, size_t> hashmap;

    for( size_t i = 0; i < vars().size(); i++ )
        hashmap[var(i).label()] = i;

    // create edge list
    vector<Edge> edges;
    edges.reserve( nrEdges );
    for( size_t i2 = 0; i2 < nrFactors(); i2++ ) {
        const VarSet& ns = factor(i2).vars();
        for( VarSet::const_iterator q = ns.begin(); q != ns.end(); q++ )
            edges.push_back( Edge(hashmap[q->label()], i2) );
    }

    // create bipartite graph
    // When dumping causal graph, sub-ids has been uniqued, so no need for check
    _G.construct( nrVars(), nrFactors(), edges.begin(), edges.end(), false );
}


/// Writes a FactorGraph to an output stream
std::ostream& dai::operator<< ( std::ostream &os, const CausalFactorGraph &fg ) {
    os << fg.nrFactors() << endl;

    for( size_t I = 0; I < fg.nrFactors(); I++ ) {
        os << endl;
        os << fg.factor(I).head().label() << endl;
        os << static_cast<char>(fg.factor(I).type) << endl;
        if (fg.factor(I).type == CausalFactor::Singleton) {
            os << setw(os.precision()+4) << fg.factor(I).prob() << endl;
        } else {
            os << fg.factor(I).body().size() << endl;
            for (const auto & i : fg.factor(I).body())
                os << i.label() << " ";
            os << endl;
        }
        // Note: node states always binary
//        for( VarSet::const_iterator i = fg.factor(I).vars().begin(); i != fg.factor(I).vars().end(); i++ )
//            os << i->states() << " ";
//        os << endl;
//        size_t nr_nonzeros = 0;
//        for( size_t k = 0; k < fg.factor(I).nrStates(); k++ )
//            if( fg.factor(I)[k] != (Real)0 )
//                nr_nonzeros++;
//        os << nr_nonzeros << endl;
//        for( size_t k = 0; k < fg.factor(I).nrStates(); k++ )
//            if( fg.factor(I)[k] != (Real)0 )
//                os << k << " " << setw(os.precision()+4) << fg.factor(I)[k] << endl;
    }

    return(os);
}


/// Reads a FactorGraph from an input stream
std::istream& dai::operator>> ( std::istream& is, CausalFactorGraph& fg ) {
    long verbose = 0;

    vector<CausalFactor> facs;
    size_t nr_Factors;
    string line;
    
    getline(is, line);
    while( line.front() == '#' )
        getline(is,line);
    nr_Factors = std::stoi(line);
    if( is.fail() )
        DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of causal factors");
    if( verbose >= 1 )
        cerr << "Reading " << nr_Factors << " causal factors..." << endl;

    for( size_t I = 0; I < nr_Factors; I++ ) {
        getline (is,line);
        if( is.fail() || line.size() > 0 )
            DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Expecting empty line");
        
        if( verbose >= 2 )
            cerr << "Reading causal factor " << I << "..." << endl;

        size_t head_id;
        getline(is, line);
        while( line.front() == '#' )
            getline(is,line);
        head_id = std::stoi(line);
        if( verbose >= 2 )
            cerr << " head: " << head_id << endl;
        Var head{head_id, 2};
        
        char type_ch;
        CausalFactor::CausalType type;
        getline(is, line);
        while( line.front() == '#' )
            getline(is,line);
        type_ch = line.front();
        switch (type_ch) {
            case CausalFactor::Singleton: {
                type = CausalFactor::Singleton;
                if( verbose >= 2 )
                    cerr << " type: singleton" << endl;
                break;
            }
            case CausalFactor::DefiniteAnd: {
                type = CausalFactor::DefiniteAnd;
                if( verbose >= 2 )
                    cerr << " type: definitive-and" << endl;
                break;
            }
            case CausalFactor::DefiniteOr: {
                type = CausalFactor::DefiniteOr;
                if( verbose >= 2 )
                    cerr << " type: definitive-or" << endl;
                break;
            }
            default: {
                DAI_THROW(UNKNOWN_ENUM_VALUE);
            }
        }

        if (type == CausalFactor::Singleton) {
            Real prob = 0;
            getline(is, line);
            while( line.front() == '#' )
                getline(is,line);
            prob = static_cast<Real>(std::stold(line));
            if( verbose >= 2 )
                cerr << " probability: " << setw(cerr.precision()+4) << prob << endl;
            facs.emplace_back(head, prob);
        } else {
            Real prob = 1;
            if (line.length() > 1) {
                prob = static_cast<Real>(std::stold(line.substr(1)));
                if( verbose >= 2 )
                    cerr << " probability: " << setw(cerr.precision()+4) << prob << endl;
            } else {
                if( verbose >= 2 )
                    cerr << " definite!" << endl;
            }
            size_t body_len;    
            getline(is, line);
            while( line.front() == '#' )
                getline(is,line);
            body_len = std::stoi(line);
            if( verbose >= 2 )
                cerr << "  body len: " << body_len << endl;


            getline(is, line);
            while( line.front() == '#' )
                getline(is,line);
            istringstream body_ss(line);
            vector<size_t> body_ids;
            for( size_t bi = 0; bi < body_len; bi++ ) {
                long body_id;
                body_ss >> body_id;
                body_ids.push_back(body_id);
            }
            if( verbose >= 2 )
                cerr << "  body: " << body_ids << endl;
            
            VarSet body;
            for( size_t body_id : body_ids ) {
                Var body_i{body_id, 2};
                body.insert(body_i);
            }
            facs.emplace_back(head, body, type == CausalFactor::DefiniteAnd, prob);
        }
    }

    if( verbose >= 3 )
        cerr << "factors:" << facs << endl;

    fg = CausalFactorGraph(facs);

    return is;
}


VarSet CausalFactorGraph::Delta( size_t i ) const {
    // calculate Markov Blanket
    VarSet Del;
    bforeach( const Neighbor &I, nbV(i) ) // for all neighboring factors I of i
        bforeach( const Neighbor &j, nbF(I) ) // for all neighboring variables j of I
            Del |= var(j);

    return Del;
}


VarSet CausalFactorGraph::Delta( const VarSet &ns ) const {
    VarSet result;
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
        result |= Delta( findVar(*n) );
    return result;
}


SmallSet<size_t> CausalFactorGraph::Deltai( size_t i ) const {
    // calculate Markov Blanket
    SmallSet<size_t> Del;
    bforeach( const Neighbor &I, nbV(i) ) // for all neighboring factors I of i
        bforeach( const Neighbor &j, nbF(I) ) // for all neighboring variables j of I
            Del |= j;

    return Del;
}


void CausalFactorGraph::makeCavity( size_t i, bool backup ) {
    // fills all Factors that include var(i) with ones
//    map<size_t,CausalFactor> newFacs;
//    bforeach( const Neighbor &I, nbV(i) ) // for all neighboring factors I of i
//        newFacs[I] = CausalFactor( factor(I).vars(), (Real)1 );
//    setFactors( newFacs, backup );
    DAI_THROW(NOT_IMPLEMENTED);
}


void CausalFactorGraph::makeRegionCavity( std::vector<size_t> facInds, bool backup ) {
//    map<size_t,Factor> newFacs;
//    for( size_t I = 0; I < facInds.size(); I++ )
//        newFacs[facInds[I]] = Factor(factor(facInds[I]).vars(), (Real)1);
//    setFactors( newFacs, backup );
    DAI_THROW(NOT_IMPLEMENTED);
}


void CausalFactorGraph::ReadFromFile( const char *filename ) {
    ifstream infile;
    infile.open( filename );
    if( infile.is_open() ) {
        infile >> *this;
        infile.close();
    } else
        DAI_THROWE(CANNOT_READ_FILE,"Cannot read from file " + std::string(filename));
}


void CausalFactorGraph::WriteToFile( const char *filename, size_t precision ) const {
    ofstream outfile;
    outfile.open( filename );
    if( outfile.is_open() ) {
        outfile.precision( precision );
        outfile << *this;
        outfile.close();
    } else
        DAI_THROWE(CANNOT_WRITE_FILE,"Cannot write to file " + std::string(filename));
}


void CausalFactorGraph::printDot( std::ostream &os ) const {
    os << "graph FactorGraph {" << endl;
    os << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
    for( size_t i = 0; i < nrVars(); i++ )
        os << "\tv" << var(i).label() << ";" << endl;
    os << "node[shape=box,width=0.3,height=0.3,fixedsize=true];" << endl;
    for( size_t I = 0; I < nrFactors(); I++ )
        os << "\tf" << I << ";" << endl;
    for( size_t i = 0; i < nrVars(); i++ )
        bforeach( const Neighbor &I, nbV(i) )  // for all neighboring factors I of i
            os << "\tv" << var(i).label() << " -- f" << I << ";" << endl;
    os << "}" << endl;
}


GraphAL CausalFactorGraph::MarkovGraph() const {
    GraphAL G( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        bforeach( const Neighbor &I, nbV(i) )
            bforeach( const Neighbor &j, nbF(I) )
                if( i < j )
                    G.addEdge( i, j, true );
    return G;
}


bool CausalFactorGraph::isMaximal( size_t I ) const {
    const VarSet& I_vars = factor(I).vars();
    size_t I_size = I_vars.size();

    if( I_size == 0 ) {
        for( size_t J = 0; J < nrFactors(); J++ )
            if( J != I )
                if( factor(J).vars().size() > 0 )
                    return false;
        return true;
    } else {
        bforeach( const Neighbor& i, nbF(I) ) {
            bforeach( const Neighbor& J, nbV(i) ) {
                if( J != I )
                    if( (factor(J).vars() >> I_vars) && (factor(J).vars().size() != I_size) )
                        return false;
            }
        }
        return true;
    }
}


size_t CausalFactorGraph::maximalFactor( size_t I ) const {
    const VarSet& I_vars = factor(I).vars();
    size_t I_size = I_vars.size();

    if( I_size == 0 ) {
        for( size_t J = 0; J < nrFactors(); J++ )
            if( J != I )
                if( factor(J).vars().size() > 0 )
                    return maximalFactor( J );
        return I;
    } else {
        bforeach( const Neighbor& i, nbF(I) ) {
            bforeach( const Neighbor& J, nbV(i) ) {
                if( J != I )
                    if( (factor(J).vars() >> I_vars) && (factor(J).vars().size() != I_size) )
                        return maximalFactor( J );
            }
        }
        return I;
    }
}


vector<VarSet> CausalFactorGraph::maximalFactorDomains() const {
    vector<VarSet> result;

    for( size_t I = 0; I < nrFactors(); I++ )
        if( isMaximal( I ) )
            result.push_back( factor(I).vars() );

    if( result.size() == 0 )
        result.push_back( VarSet() );
    return result;
}


Real CausalFactorGraph::logScore( const std::vector<size_t>& statevec ) const {
    // Construct a State object that represents statevec
    // This decouples the representation of the joint state in statevec from the factor graph
    map<Var, size_t> statemap;
    for( size_t i = 0; i < statevec.size(); i++ )
        statemap[var(i)] = statevec[i];
    State S(statemap);

    // Evaluate the log probability of the joint configuration in statevec
    // by summing the log factor entries of the factors that correspond to this joint configuration
    Real lS = 0.0;
    for( size_t I = 0; I < nrFactors(); I++ )
        if( factor(I).type == CausalFactor::Singleton ) {
            lS += dai::log( 1-factor(I).prob() );
            lS += dai::log( factor(I).prob() );
        }
    return lS;
}


void CausalFactorGraph::clamp( size_t i, size_t x, bool backup ) {
    DAI_ASSERT( x <= var(i).states() );

    map<size_t, CausalFactor> newFacs;
    bforeach( const Neighbor &I, nbV(i) ) {
        CausalFactor newFac = factor(I).gen_clamped(var(i), x);
        newFacs.emplace(I, newFac);
    }
    setFactors( newFacs, backup );

    return;
}


void CausalFactorGraph::clampVar( size_t i, const vector<size_t> &is, bool backup ) {
    if (is.size() == 1)
        clamp(i, is[0], backup);
    else
        DAI_THROW(NOT_IMPLEMENTED);
}


void CausalFactorGraph::clampFactor( size_t I, const vector<size_t> &is, bool backup ) {
    DAI_THROW(NOT_IMPLEMENTED);
}


void CausalFactorGraph::backupFactor( size_t I ) {
    map<size_t,CausalFactor>::iterator it = _backup.find( I );
    if( it != _backup.end() )
        DAI_THROW(MULTIPLE_UNDO);
    _backup.emplace(I, factor(I));
}


void CausalFactorGraph::restoreFactor( size_t I ) {
    map<size_t,CausalFactor>::iterator it = _backup.find( I );
    if( it != _backup.end() ) {
        setFactor(I, it->second);
        _backup.erase(it);
    } else
        DAI_THROW(OBJECT_NOT_FOUND);
}


void CausalFactorGraph::backupFactors( const VarSet &ns ) {
    for( size_t I = 0; I < nrFactors(); I++ )
        if( factor(I).vars().intersects( ns ) )
            backupFactor( I );
}


void CausalFactorGraph::restoreFactors( const VarSet &ns ) {
    map<size_t,CausalFactor> facs;
    for( map<size_t,CausalFactor>::iterator uI = _backup.begin(); uI != _backup.end(); ) {
        if( factor(uI->first).vars().intersects( ns ) ) {
            facs.insert( *uI );
            _backup.erase(uI++);
        } else
            uI++;
    }
    setFactors( facs );
}


void CausalFactorGraph::restoreFactors() {
    setFactors( _backup );
    _backup.clear();
}


void CausalFactorGraph::backupFactors( const std::set<size_t> & facs ) {
    for( std::set<size_t>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ )
        backupFactor( *fac );
}


bool CausalFactorGraph::isPairwise() const {
    bool pairwise = true;
    for( size_t I = 0; I < nrFactors() && pairwise; I++ )
        if( factor(I).vars().size() > 2 )
            pairwise = false;
    return pairwise;
}

CausalFactorGraph CausalFactorGraph::clamped( size_t i, size_t state ) const {
    DAI_THROW(NOT_IMPLEMENTED);
    return CausalFactorGraph( *this );
}


CausalFactorGraph CausalFactorGraph::maximalFactors() const {
    DAI_THROW(NOT_IMPLEMENTED);
    return CausalFactorGraph( *this );
}
