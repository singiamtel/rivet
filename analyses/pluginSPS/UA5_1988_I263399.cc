// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  /// @brief UA5 charged particle correlations at 200, 546 and 900 GeV
  class UA5_1988_I263399 : public Analysis {
  public:

    UA5_1988_I263399()
      : Analysis("UA5_1988_I263399")
    {    }


    /// @name Analysis methods
    /// @{

    void init() {
      // Projections
      declare(TriggerUA5(), "Trigger");

      // Symmetric eta interval
      declare(ChargedFinalState((Cuts::etaIn(-0.5, 0.5))), "CFS05");

      // Asymmetric intervals first
      // Forward eta intervals
      declare(ChargedFinalState((Cuts::etaIn(0.0, 1.0))), "CFS10F");
      declare(ChargedFinalState((Cuts::etaIn(0.5, 1.5))), "CFS15F");
      declare(ChargedFinalState((Cuts::etaIn(1.0, 2.0))), "CFS20F");
      declare(ChargedFinalState((Cuts::etaIn(1.5, 2.5))), "CFS25F");
      declare(ChargedFinalState((Cuts::etaIn(2.0, 3.0))), "CFS30F");
      declare(ChargedFinalState((Cuts::etaIn(2.5, 3.5))), "CFS35F");
      declare(ChargedFinalState((Cuts::etaIn(3.0, 4.0))), "CFS40F");

      // Backward eta intervals
      declare(ChargedFinalState((Cuts::etaIn(-1.0,  0.0))), "CFS10B");
      declare(ChargedFinalState((Cuts::etaIn(-1.5, -0.5))), "CFS15B");
      declare(ChargedFinalState((Cuts::etaIn(-2.0, -1.0))), "CFS20B");
      declare(ChargedFinalState((Cuts::etaIn(-2.5, -1.5))), "CFS25B");
      declare(ChargedFinalState((Cuts::etaIn(-3.0, -2.0))), "CFS30B");
      declare(ChargedFinalState((Cuts::etaIn(-3.5, -2.5))), "CFS35B");
      declare(ChargedFinalState((Cuts::etaIn(-4.0, -3.0))), "CFS40B");

      // Histogram booking, we have sqrt(s) = 200, 546 and 900 GeV
      // TODO use Estimate1D to be able to output errors
      if (isCompatibleWithSqrtS(200.0*GeV)) {
        book(_hist_correl, 2, 1, 1);
        book(_hist_correl_asym, 3, 1, 1);
      } else if (isCompatibleWithSqrtS(546.0*GeV)) {
        book(_hist_correl, 2, 1, 2);
        book(_hist_correl_asym, 3, 1, 2);
      } else if (isCompatibleWithSqrtS(900.0*GeV)) {
        book(_hist_correl, 2, 1, 3);
        book(_hist_correl_asym, 3, 1, 3);
      }

      book(_sumWPassed, "sumW");
    }


    void analyze(const Event& event) {
      // Trigger
      const bool trigger = apply<TriggerUA5>(event, "Trigger").nsdDecision();
      if (!trigger) vetoEvent;
      _sumWPassed->fill();

      // Count forward/backward particles
      n_10f.push_back(apply<ChargedFinalState>(event, "CFS10F").size());
      n_15f.push_back(apply<ChargedFinalState>(event, "CFS15F").size());
      n_20f.push_back(apply<ChargedFinalState>(event, "CFS20F").size());
      n_25f.push_back(apply<ChargedFinalState>(event, "CFS25F").size());
      n_30f.push_back(apply<ChargedFinalState>(event, "CFS30F").size());
      n_35f.push_back(apply<ChargedFinalState>(event, "CFS35F").size());
      n_40f.push_back(apply<ChargedFinalState>(event, "CFS40F").size());
      //
      n_10b.push_back(apply<ChargedFinalState>(event, "CFS10B").size());
      n_15b.push_back(apply<ChargedFinalState>(event, "CFS15B").size());
      n_20b.push_back(apply<ChargedFinalState>(event, "CFS20B").size());
      n_25b.push_back(apply<ChargedFinalState>(event, "CFS25B").size());
      n_30b.push_back(apply<ChargedFinalState>(event, "CFS30B").size());
      n_35b.push_back(apply<ChargedFinalState>(event, "CFS35B").size());
      n_40b.push_back(apply<ChargedFinalState>(event, "CFS40B").size());
      //
      n_05 .push_back(apply<ChargedFinalState>(event, "CFS05").size());
    }


    void finalize() {
      // The correlation strength is defined in formulas
      // 4.1 and 4.2

      // Fill histos, gap width histo comes first
      //      * Set the errors as Delta b / sqrt(sumWPassed) with
      //      Delta b being the absolute uncertainty of b according to
      //      Gaussian error-propagation (linear limit) and assuming
      //      Poissonian uncertainties for the number of particles in
      //      the eta-intervals
      //

      // Fill the DPS
      pair<double,double> corr;
      corr = correlation_helper(n_10f, n_10b, _sumWPassed);
      _hist_correl->bin(1).set(corr.first, corr.second);
      corr = correlation_helper(n_15f, n_15b, _sumWPassed);
      _hist_correl->bin(2).set(corr.first, corr.second);
      corr = correlation_helper(n_20f, n_20b, _sumWPassed);
      _hist_correl->bin(3).set(corr.first, corr.second);
      corr = correlation_helper(n_25f, n_25b, _sumWPassed);
      _hist_correl->bin(4).set(corr.first, corr.second);
      corr = correlation_helper(n_30f, n_30b, _sumWPassed);
      _hist_correl->bin(5).set(corr.first, corr.second);
      corr = correlation_helper(n_35f, n_35b, _sumWPassed);
      _hist_correl->bin(6).set(corr.first, corr.second);
      corr = correlation_helper(n_40f, n_40b, _sumWPassed);
      _hist_correl->bin(7).set(corr.first, corr.second);

      // Fill gap-center histo (Fig 15)
      //
      // The first bin contains the c_str strengths of
      // the gap size histo that has ane eta gap of two
      //
      // Fill in correlation strength for assymetric intervals,
      // see Tab. 5
      // Fill the DPS
      corr = correlation_helper(n_20f, n_20b, _sumWPassed);
      _hist_correl_asym->bin(1).set(corr.first, corr.second);
      corr = correlation_helper(n_25f, n_15b, _sumWPassed);
      _hist_correl_asym->bin(2).set(corr.first, corr.second);
      corr = correlation_helper(n_30f, n_10b, _sumWPassed);
      _hist_correl_asym->bin(3).set(corr.first, corr.second);
      corr = correlation_helper(n_35f, n_05 , _sumWPassed);
      _hist_correl_asym->bin(4).set(corr.first, corr.second);
      corr = correlation_helper(n_40f, n_10f, _sumWPassed);
      _hist_correl_asym->bin(5).set(corr.first, corr.second);

    }

    /// @}


  private:

    /// Helper function to fill correlation points into scatter plot
    pair<double,double> correlation_helper(const vector<int>& nf, const vector<int>& nb, CounterPtr sumWPassed) {
      return make_pair(correlation(nf, nb),  correlation_err(nf, nb)/sqrt(sumWPassed->val()));
    }

    /// Counter
    CounterPtr _sumWPassed;


    /// @name Vectors for storing the number of particles in the different eta intervals per event.
    /// @todo A better way is needed to make this re-entrant
    /// @{
    vector<int> n_10f, n_15f, n_20f, n_25f, n_30f, n_35f, n_40f;
    vector<int> n_10b, n_15b, n_20b, n_25b, n_30b, n_35b, n_40b;
    vector<int> n_05;
    /// @}


    /// @name Histograms
    /// @{
    // Symmetric eta intervals
    Estimate1DPtr _hist_correl;
    // For asymmetric eta intervals
    Estimate1DPtr _hist_correl_asym;
    /// @}

  };


  RIVET_DECLARE_ALIASED_PLUGIN(UA5_1988_I263399, UA5_1988_S1867512);

}
