//
// Created by Yifan Chen on 2023/8/21.
//

#ifndef LIBDAI_MULTITHREAD_EMALG_H
#define LIBDAI_MULTITHREAD_EMALG_H

#include "dai/alldai.h"

namespace dai {

    /// MultiEMAlg performs a multi-threaded version of EM, which performs E-steps in parallel for each evidence.
    /** Need to specify max_threads, plus Evidence / InfAlg / MSteps
     *
     *  \author Yifan Chen
     */
    class MultiEMAlg {
    private:
        /// All the data samples used during learning
        const Evidence &_evidence;

        /// How to do the expectation step
        InfAlg &_estep;

        /// The maximization steps to take
        std::vector<MaximizationStep> _msteps;

        /// Number of iterations done
        size_t _iters;

        /// History of likelihoods
        std::vector<Real> _lastLogZ;

        /// Maximum number of iterations
        size_t _max_iters;

        /// Convergence tolerance
        Real _log_z_tol;

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
        /// Maximum number of threads for paralle updating
        static const size_t MAX_THREADS;

        /// Construct an EMAlg from several objects
        /** \param evidence Specifies the observed evidence
         *  \param estep Inference algorithm to be used for the E-step
         *  \param msteps Vector of maximization steps, each of which is a group of parameter estimation tasks
         *  \param termconditions Termination conditions @see setTermConditions()
         *  \param threads max threads for parallel estep
         */
        MultiEMAlg(const Evidence &evidence, InfAlg &estep, std::vector<MaximizationStep> &msteps,
                   const PropertySet &termconditions)
                : _evidence(evidence), _estep(estep), _msteps(msteps), _iters(0), _lastLogZ(),
                  _max_iters(MAX_ITERS_DEFAULT), _log_z_tol(LOG_Z_TOL_DEFAULT), _max_jobs(MAX_THREADS) {
            setTermConditions(termconditions);
        }

        /// Construct an EMAlg from Evidence \a evidence, an InfAlg \a estep, and an input stream \a mstep_file
        /** \see \ref fileformats-emalg
         */
        MultiEMAlg(const Evidence &evidence, InfAlg &estep, std::istream &mstep_file);

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
        const InfAlg &eStep() const { return _estep; }

        /// Iterate once over all maximization steps
        /** \return Log-likelihood after iteration
         */
        Real iterate();

        /// Iterate over a single MaximizationStep
        Real iterate(MaximizationStep &mstep);

        /// Iterate until termination conditions are satisfied
        void run();

        /// \name Iterator interface
        //@{
        /// Iterator over the maximization steps
        typedef std::vector<MaximizationStep>::iterator s_iterator;
        /// Constant iterator over the maximization steps
        typedef std::vector<MaximizationStep>::const_iterator const_s_iterator;

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
}
#endif //LIBDAI_MULTITHREAD_EMALG_H
