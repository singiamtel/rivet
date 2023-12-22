// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief ALEPH LEP1 charged multiplicity in hadronic Z decay
  ///
  /// @author Andy Buckley
  class ALEPH_1991_I319520 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_1991_I319520);


    /// @name Analysis methods
    /// @{

    /// Book projections and histogram
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      book(_histChTot ,1, 1, 1);
      book(_histAver, 2, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histChTot->fill(cfs.size());
      _histAver->fill(Ecms, cfs.size());
    }


    /// Normalize the histogram
    void finalize() {
      scale(_histChTot, 1./sumOfWeights()); // same as in ALEPH 1996
      scale(_histAver , 1./sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _histChTot;
    BinnedHistoPtr<string> _histAver;
    const string Ecms = "91.25";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ALEPH_1991_I319520, ALEPH_1991_S2435284);

}
