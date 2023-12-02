// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 Z + jet + \f$ X \f$ cross-section / \f$ p_\perp \f$ distributions
  class D0_2009_S8202443 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2009_S8202443);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      // Leptons in constrained tracking acceptance
      Cut cuts = (Cuts::abseta < 1.1 || Cuts::absetaIn(1.5, 2.5)) && Cuts::pT > 25*GeV;
      DileptonFinder zfinder_constrained(91.2*GeV, 0.2, cuts &&
                                         Cuts::abspid == PID::ELECTRON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zfinder_constrained, "DileptonFinderConstrained");
      FastJets conefinder_constrained(zfinder_constrained.remainingFinalState(), JetAlg::D0ILCONE, 0.5);
      declare(conefinder_constrained, "ConeFinderConstrained");

      // Unconstrained leptons
      DileptonFinder zfinder(91.2*GeV, 0.2, Cuts::abspid == PID::ELECTRON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zfinder, "DileptonFinder");
      FastJets conefinder(zfinder.remainingFinalState(), JetAlg::D0ILCONE, 0.5);
      declare(conefinder, "ConeFinder");

      book(_h_jet1_pT_constrained ,1, 1, 1);
      book(_h_jet2_pT_constrained ,3, 1, 1);
      book(_h_jet3_pT_constrained ,5, 1, 1);
      book(_h_jet1_pT ,2, 1, 1);
      book(_h_jet2_pT ,4, 1, 1);
      book(_h_jet3_pT ,6, 1, 1);

      book(_sum_of_weights,"sum_of_weights");
      book(_sum_of_weights_constrained, "sum_of_weights_constrained");
    }


    // Do the analysis
    void analyze(const Event& e) {
      // Unconstrained electrons
      const DileptonFinder& zfinder = apply<DileptonFinder>(e, "DileptonFinder");
      if (zfinder.bosons().size() == 0) {
        MSG_DEBUG("No unique lepton pair found.");
        vetoEvent;
      }
      _sum_of_weights->fill();
      const Jets jets_cut = apply<JetFinder>(e, "ConeFinder").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.5);
      if (jets_cut.size() > 0)
        _h_jet1_pT->fill(jets_cut[0].pT()/GeV);
      if (jets_cut.size() > 1)
        _h_jet2_pT->fill(jets_cut[1].pT()/GeV);
      if (jets_cut.size() > 2)
        _h_jet3_pT->fill(jets_cut[2].pT()/GeV);


      // Constrained electrons
      const DileptonFinder& zfinder_constrained = apply<DileptonFinder>(e, "DileptonFinderConstrained");
      if (zfinder_constrained.bosons().size() == 0) {
        MSG_DEBUG("No unique constrained lepton pair found.");
        return; // Not really a "veto", since if we got this far there is an unconstrained Z
      }
      _sum_of_weights_constrained->fill();
      const Jets& jets_constrained = apply<JetFinder>(e, "ConeFinderConstrained").jetsByPt(Cuts::pT > 20*GeV);
      /// @todo Replace this explicit selection with a Cut
      Jets jets_cut_constrained;
      for (const Jet& j : jets_constrained) {
        if (j.abseta() < 2.5) jets_cut_constrained.push_back(j);
      }
      if (jets_cut_constrained.size() > 0)
        _h_jet1_pT_constrained->fill(jets_cut_constrained[0].pT()/GeV);
      if (jets_cut_constrained.size() > 1)
        _h_jet2_pT_constrained->fill(jets_cut_constrained[1].pT()/GeV);
      if (jets_cut_constrained.size() > 2)
        _h_jet3_pT_constrained->fill(jets_cut_constrained[2].pT()/GeV);
    }


    // Finalize
    void finalize() {
      scale(_h_jet1_pT, 1/ *_sum_of_weights);
      scale(_h_jet2_pT, 1/ *_sum_of_weights);
      scale(_h_jet3_pT, 1/ *_sum_of_weights);
      scale(_h_jet1_pT_constrained, 1/ *_sum_of_weights_constrained);
      scale(_h_jet2_pT_constrained, 1/ *_sum_of_weights_constrained);
      scale(_h_jet3_pT_constrained, 1/ *_sum_of_weights_constrained);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_jet1_pT;
    Histo1DPtr _h_jet2_pT;
    Histo1DPtr _h_jet3_pT;
    Histo1DPtr _h_jet1_pT_constrained;
    Histo1DPtr _h_jet2_pT_constrained;
    Histo1DPtr _h_jet3_pT_constrained;
    /// @}

    CounterPtr _sum_of_weights, _sum_of_weights_constrained;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2009_S8202443, D0_2009_I815094);

}
