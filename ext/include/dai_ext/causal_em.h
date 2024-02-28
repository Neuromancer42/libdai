//
// Created by Yifan Chen on 2024/1/14.
//

#ifndef LIBDAI_CAUSAL_EM_H
#define LIBDAI_CAUSAL_EM_H

#include "dai/alldai.h"
#include "dai_ext/causal_bp.h"
#include "dai_ext/causal_fg.h"

namespace dai {
    
    /// Represents a single factor or set of factors whose parameters should be estimated.
    /** To ensure that parameters can be shared between different factors during
     *  EM learning, each factor's values are reordered to match a desired variable
     *  ordering. The ordering of the variables in a factor may therefore differ
     *  from the canonical ordering used in libDAI. The SharedParameters
     *  class combines one or more factors (together with the specified orderings
     *  of the variables) with a ParameterEstimation object, taking care of the
     *  necessary permutations of the factor entries / parameters.
     * 
     *  \author Charles Vaske
     */
    class CausalSharedParam {
    public:
        /// Convenience label for an index of a factor in a FactorGraph.
        typedef size_t FactorIndex;
        /// Convenience label for a grouping of factor orientations.
        typedef std::map<FactorIndex, Var > FactorOrientations;

    private:
        /// Maps factor indices to the corresponding desired variable orderings
        FactorOrientations _varorders;
        /// Parameter estimation method to be used
        ParameterEstimation *_estimation;
        /// Indicates whether \c *this gets ownership of _estimation
        bool _ownEstimation;
        /// The accumulated expectations
        Prob* _expectations;
        
        /// Initializes _varsets and _perms from _varorders and checks whether their state spaces correspond with _estimation.probSize()
        void setPermsAndVarSetsFromVarOrders();

    public:
        /// Constructor
        /** \param varorders  all the factor orientations for this parameter
         *  \param estimation a pointer to the parameter estimation method
         *  \param ownPE whether the constructed object gets ownership of \a estimation
         */
        CausalSharedParam( const FactorOrientations &varorders, ParameterEstimation *estimation, bool ownPE=false );

        /// Construct a SharedParameters object from an input stream \a is and a factor graph \a fg
        /** \see \ref fileformats-emalg-sharedparameters
         *  \throw INVALID_EMALG_FILE if the input stream is not valid
         */
        CausalSharedParam( std::istream &is, const CausalFactorGraph &causal_fg );

        /// Copy constructor
        CausalSharedParam( const CausalSharedParam &sp ) : _varorders(sp._varorders), _estimation(sp._estimation), _ownEstimation(sp._ownEstimation), _expectations(NULL) {
            // If sp owns its _estimation object, we should clone it instead of copying the pointer
            if( _ownEstimation )
                _estimation = _estimation->clone();
            _expectations = new Prob(*sp._expectations);
        }

        /// Destructor
        ~CausalSharedParam() {
            // If we own the _estimation object, we should delete it now
            if( _ownEstimation )
                delete _estimation;
            if( _expectations != NULL)
                delete _expectations;
        }

        /// Collect the expected values (beliefs) according to \a alg
        /** For each of the relevant factors (that shares the parameters we are interested in),
         *  the corresponding belief according to \a alg is obtained and its entries are permuted
         *  such that their ordering corresponds with the shared parameters that we are estimating.
         */
        void collectExpectations( CausalBP &alg );

        /// Return the current accumulated expectations
        const Prob& currentExpectations() const { return *_expectations; }

        ParameterEstimation& getPEst() const { return *_estimation; }

        /// Estimate and set the shared parameters
        /** Based on the expectation statistics collected so far, the shared parameters are estimated
         *  using the parameter estimation subclass method estimate(). Then, each of the relevant
         *  factors in \a fg (that shares the parameters we are interested in) is set according 
         *  to those parameters (permuting the parameters accordingly).
         */
        void setParameters( CausalFactorGraph &fg );

        /// Return a reference to the vector of factor orientations
        /** This is necessary for determing which variables were used
         *  to estimate parameters, and analysis of expectations
         *  after an Estimation step has been performed.
         */
        const FactorOrientations& getFactorOrientations() const { return _varorders; }

        /// Reset the current expectations
        void clear( ) { _expectations->fill(0); }
    };


/// A MaximizationStep groups together several parameter estimation tasks (SharedParameters objects) into a single unit.
/** \author Charles Vaske
 */
    class CausalMaxStep {
    private:
        /// Vector of parameter estimation tasks of which this maximization step consists
        std::vector<CausalSharedParam> _params;

    public:
        /// Default constructor
        CausalMaxStep() : _params() {}

        /// Construct MaximizationStep from a vector of parameter estimation tasks
        CausalMaxStep( std::vector<CausalSharedParam> &maximizations ) : _params(maximizations) {}

        /// Constructor from an input stream and a corresponding factor graph
        /** \see \ref fileformats-emalg-maximizationstep
         */
        CausalMaxStep( std::istream &is, const CausalFactorGraph &fg_varlookup );

        /// Collect the beliefs from this InfAlg as expectations for the next Maximization step
        void addExpectations( CausalBP &alg );

        /// Using all of the currently added expectations, make new factors with maximized parameters and set them in the FactorGraph.
        void maximize( CausalFactorGraph &fg );

        /// Clear the step, to be called at the begining of each step
        void clear( );

        /// \name Iterator interface
        //@{
        /// Iterator over the parameter estimation tasks
        typedef std::vector<CausalSharedParam>::iterator iterator;
        /// Constant iterator over the parameter estimation tasks
        typedef std::vector<CausalSharedParam>::const_iterator const_iterator;

        /// Returns iterator that points to the first parameter estimation task
        iterator begin() { return _params.begin(); }
        /// Returns constant iterator that points to the first parameter estimation task
        const_iterator begin() const { return _params.begin(); }
        /// Returns iterator that points beyond the last parameter estimation task
        iterator end() { return _params.end(); }
        /// Returns constant iterator that points beyond the last parameter estimation task
        const_iterator end() const { return _params.end(); }
        //@}
    };
    
    class CausalEM {
    private:
        /// All the data samples used during learning
        const Evidence &_evidence;

        /// How to do the expectation step
        CausalBP &_estep;

        /// The maximization steps to take
        std::vector<CausalMaxStep> _msteps;

        /// Number of iterations done
        size_t _iters;

        /// History of likelihoods
        std::vector<Real> _lastLogZ;

        /// Maximum number of iterations
        size_t _max_iters;

        /// Convergence tolerance
        Real _log_z_tol;
        
        /// Number of parallel jobs
        size_t _max_jobs;

    public:
        /// Key for setting maximum iterations
        static const std::string MAX_ITERS_KEY;
        /// Default maximum iterations
        static const size_t MAX_ITERS_DEFAULT;
        /// Key for setting likelihood termination condition
        static const std::string LOG_Z_TOL_KEY;
        /// Default likelihood tolerance
        static const Real LOG_Z_TOL_DEFAULT;
        /// Maximum number of threads for parallel updating
        static const size_t MAX_THREADS;
        
        /// Construct an EMAlg from several objects
        /** \param evidence Specifies the observed evidence
         *  \param estep Inference algorithm to be used for the E-step
         *  \param msteps Vector of maximization steps, each of which is a group of parameter estimation tasks
         *  \param termconditions Termination conditions @see setTermConditions()
         *  \param threads max threads for parallel estep
         */
        CausalEM(const Evidence &evidence, CausalBP &estep, std::vector<CausalMaxStep> &msteps,
              const PropertySet &termconditions)
                : _evidence(evidence), _estep(estep), _msteps(msteps), _iters(0), _lastLogZ(),
                  _max_iters(MAX_ITERS_DEFAULT), _log_z_tol(LOG_Z_TOL_DEFAULT), _max_jobs(MAX_THREADS) {
            setTermConditions(termconditions);
        }

        /// Construct an EMAlg from Evidence \a evidence, an CausalBP \a estep, and an input stream \a mstep_file
        /** \see \ref fileformats-emalg
         */
        CausalEM(const Evidence &evidence, CausalBP &estep, std::istream &mstep_file);

        /// Change the conditions for termination
        /** There are two possible parameters in the PropertySet \a p:
         *    - \a max_iters maximum number of iterations
         *    - \a log_z_tol critical proportion of increase in logZ
         *
         *  \see hasSatisifiedTermConditions()
         */
        void setTermConditions(const PropertySet &p);

        /// Set the number of max therads
        void setMaxJobs(size_t max_jobs) {
            std::clog << "LibDAI: set max num of jobs from " << _max_jobs << " to " << max_jobs << std::endl;
            _max_jobs = max_jobs;
        }
        
        /// Determine if the termination conditions have been met.
        /** There are two sufficient termination conditions:
         *    -# the maximum number of iterations has been performed
         *    -# the ratio of logZ increase over previous logZ is less than the
         *       tolerance, i.e.,
         *       \f$ \frac{\log(Z_t) - \log(Z_{t-1})}{| \log(Z_{t-1}) | } < \mathrm{tol} \f$.
         */
        bool hasSatisfiedTermConditions() const;

        /// Return the last calculated log likelihood
        Real logZ() const { return _lastLogZ.back(); }

        /// Returns number of iterations done so far
        size_t Iterations() const { return _iters; }

        /// Get the E-step method used
        const CausalBP &eStep() const { return _estep; }

        /// Iterate once over all maximization steps
        /** \return Log-likelihood after iteration
         */
        Real iterate();

        /// Iterate over a single MaximizationStep
        Real iterate(CausalMaxStep &mstep);

        /// Iterate until termination conditions are satisfied
        void run();

        /// \name Iterator interface
        //@{
        /// Iterator over the maximization steps
        typedef std::vector<CausalMaxStep>::iterator s_iterator;
        /// Constant iterator over the maximization steps
        typedef std::vector<CausalMaxStep>::const_iterator const_s_iterator;

        /// Returns iterator that points to the first maximization step
        s_iterator s_begin() { return _msteps.begin(); }

        /// Returns constant iterator that points to the first maximization step
        const_s_iterator s_begin() const { return _msteps.begin(); }

        /// Returns iterator that points beyond the last maximization step
        s_iterator s_end() { return _msteps.end(); }

        /// Returns constant iterator that points beyond the last maximization step
        const_s_iterator s_end() const { return _msteps.end(); }
        //@}
    };
} // end of namespace dai

#endif //LIBDAI_CAUSAL_EM_H
