// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Average multiplcity at a range of energies
  class JADE_1983_I190818 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1983_I190818);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      if( !(isCompatibleWithSqrtS(12.0*GeV) ||
	    isCompatibleWithSqrtS(30.0*GeV) ||
	    isCompatibleWithSqrtS(35.0*GeV) )) {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      book(_counter, "/TMP/MULT");
      book(_mult, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _counter->fill(cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_counter,1./sumOfWeights());

      double val = _counter->val();
      double err = _counter->err();

      for (auto& b : _mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(val, err);
        }
      }
    }
    /// @}

  private:

    // Histogram
    CounterPtr _counter;
    Estimate1DPtr _mult;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(JADE_1983_I190818);

}
