// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2011_S8971293 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_S8971293);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      /// Initialise and register projections
      declare(FastJets(FinalState(), JetAlg::ANTIKT, 0.6), "AntiKtJets06");

      /// Book histograms
      book(_h_deltaPhi, {110., 160., 210., 260., 310., 400., 500., 600., 800., 10000.});
      for (auto& b : _h_deltaPhi->bins()) {
        book(b, 1, 1, b.index());
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Jets jets06 = apply<FastJets>(event, "AntiKtJets06").jetsByPt(Cuts::absrap < 2.8 && Cuts::pT > 100*GeV);
      if (jets06.size()>1){
        if (jets06[0].absrap()<0.8 && jets06[1].absrap()<0.8) {
          double observable = mapAngle0ToPi(jets06[0].phi()-jets06[1].phi()) / M_PI;
          _h_deltaPhi->fill(jets06[0].pT(), observable);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_deltaPhi, 1.0/M_PI);
    }

    /// @}


  private:

    /// Histograms
    Histo1DGroupPtr _h_deltaPhi;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_S8971293, ATLAS_2011_I889546);

}
