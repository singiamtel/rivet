// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// Rivet analysis class for ATLAS_2015_I1387176 dataset
  class ATLAS_2015_I1387176 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2015_I1387176);


    /// Initialization, called once before running
    void init() {
      // Projections
      FastJets jets(FinalState(), FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "Jets");

      // Book histograms
      book(_hist_EEC,  1, 1, 1);
      book(_hist_AEEC, 2, 1, 1);

      // add dummy histogram for heterogenous merging
      // @todo What is this for exactly? o.0
      book(_hist_dummy, "d01-x01-y02", refData("d01-x01-y01"));
    }

    void analyze(const Event& event) {

      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 50.0*GeV && Cuts::abseta < 2.5);

      if (jets.size() < 2)  vetoEvent;
      if (jets[0].pT() + jets[1].pT() < 500*GeV)  vetoEvent;

      double sumEt = 0.0;
      for (Jet j : jets)  sumEt += j.E() / cosh(j.eta());

      for (Jet j1 : jets) {
        double et1 = j1.E() / cosh(j1.eta());

        for (Jet j2 : jets) {
          double et2 = j2.E() / cosh(j2.eta());
          double etWeight = et1 * et2 / ( sumEt * sumEt );
          double dPhi = deltaPhi(j1, j2);
          double cosPhi = cos(dPhi);
          if (cosPhi == 1.0)  cosPhi = 0.9999;

          _hist_EEC->fill(cosPhi, etWeight);
          _hist_dummy->fill(cosPhi, etWeight);
	      }
      }
    }

    void finalize() {

      scale(_hist_dummy, crossSectionPerEvent());
      normalize(_hist_EEC);

      size_t nBins = _hist_EEC->numBins();
      for (size_t k = 1; k < (nBins/2)+1; ++k) {
        const double dV = _hist_EEC->bin(k).dVol();
        const double y = (_hist_EEC->bin(k).sumW() - _hist_EEC->bin(nBins-k+1).sumW())/dV;
        const double e1 = _hist_EEC->bin(k).errW()/dV;
        const double e2 = _hist_EEC->bin(nBins-k+1).errW()/dV;
        const double ey = sqrt( e1 * e1 + e2 * e2 );
        _hist_AEEC->bin(k).set(y, ey);
      }

    }

  private:
    Histo1DPtr _hist_EEC;
    Histo1DPtr _hist_dummy;
    Estimate1DPtr _hist_AEEC;
  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2015_I1387176);

}
