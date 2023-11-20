// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {



  /// @brief MC validation analysis for jet events
  class MC_JETTAGS : public Analysis {
  public:

    MC_JETTAGS()
      : Analysis("MC_JETTAGS")
    {    }


    void init() {
      FinalState fs;
      declare(FastJets(fs, JetAlg::ANTIKT, 0.4), "Jets04");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.6), "Jets06");

      book(_h_numBTagsPerJet[0] ,"numBTagsPer04Jet", 5, -0.5, 4.5);
      book(_h_numBTagsPerJet[1] ,"numBTagsPer06Jet", 5, -0.5, 4.5);
      book(_h_numCTagsPerJet[0] ,"numCTagsPer04Jet", 5, -0.5, 4.5);
      book(_h_numCTagsPerJet[1] ,"numCTagsPer06Jet", 5, -0.5, 4.5);
      book(_h_numTauTagsPerJet[0] ,"numTauTagsPer04Jet", 5, -0.5, 4.5);
      book(_h_numTauTagsPerJet[1] ,"numTauTagsPer06Jet", 5, -0.5, 4.5);
    }


    void analyze(const Event& event) {

      const Jets jets04 = apply<FastJets>(event, "Jets04").jetsByPt(Cuts::pT > 20*GeV);
      const Jets jets06 = apply<FastJets>(event, "Jets06").jetsByPt(Cuts::pT > 20*GeV);

      for (const Jet& j : jets04) {
        _h_numBTagsPerJet[0]->fill(j.bTags().size());
        _h_numCTagsPerJet[0]->fill(j.cTags().size());
        _h_numTauTagsPerJet[0]->fill(j.tauTags().size());
      }
      for (const Jet& j : jets06) {
        _h_numBTagsPerJet[1]->fill(j.bTags().size());
        _h_numCTagsPerJet[1]->fill(j.cTags().size());
        _h_numTauTagsPerJet[1]->fill(j.tauTags().size());
      }
    }


    void finalize() {
      for (size_t i = 0; i < 2; ++i) {
        normalize(_h_numBTagsPerJet[i]);
        normalize(_h_numCTagsPerJet[i]);
        normalize(_h_numTauTagsPerJet[i]);
      }
    }


  private:

    Histo1DPtr _h_numBTagsPerJet[2], _h_numCTagsPerJet[2], _h_numTauTagsPerJet[2];

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_JETTAGS);

}
