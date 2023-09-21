// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// CMS azimuthal decorrelations
  class CMS_2011_S8950903 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8950903);


    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.5);
      declare(akt, "antikT");

      book(_h_deltaPhi, {80., 110., 140., 200., 300., 7000.},
                        {"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01", "d05-x01-y01"});
    }


    void analyze(const Event & event) {

      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      if (fabs(jets[0].eta()) > 1.1 || jets[0].pT() < 80.) vetoEvent;
      if (fabs(jets[1].eta()) > 1.1 || jets[1].pT() < 30.) vetoEvent;

      double dphi = deltaPhi(jets[0].momentum(), jets[1].phi());

      _h_deltaPhi->fill(jets[0].pT(), dphi);
    }


    void finalize() {
      normalize(_h_deltaPhi);
    }


  private:

    Histo1DGroupPtr _h_deltaPhi;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8950903, CMS_2011_I885663);

}
