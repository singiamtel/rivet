// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Multijet transverse energy-energy correlations (TEEC) at 8 TeV
  class ATLAS_2017_I1609253 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2017_I1609253);


    /// Initialization, called once before running
    void init() {

      // Projections
      const FastJets jets(FinalState(), JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::ALL);
      declare(jets, "Jets");

      // Book histograms
      book(_hist_EEC1 ,  1, 1, 1);
      book(_hist_AEEC1,  2, 1, 1);
      book(_hist_EEC2 ,  3, 1, 1);
      book(_hist_AEEC2,  4, 1, 1);
      book(_hist_EEC3 ,  5, 1, 1);
      book(_hist_AEEC3,  6, 1, 1);
      book(_hist_EEC4 ,  7, 1, 1);
      book(_hist_AEEC4,  8, 1, 1);
      book(_hist_EEC5 ,  9, 1, 1);
      book(_hist_AEEC5, 10, 1, 1);
      book(_hist_EEC6 , 11, 1, 1);
      book(_hist_AEEC6, 12, 1, 1);
    }


    void analyze(const Event& event) {

      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 100*GeV);
      if (jets.size() < 2)  vetoEvent;

      double sumPt12 = jets[0].pt() + jets[1].pt();
      if (sumPt12 < 800*GeV)  vetoEvent;

      double sumEt = 0.;
      for (const Jet& j : jets) sumEt += j.Et();

      for (const Jet& j1 : jets) {
        double et1 = j1.Et();

        for (const Jet& j2 : jets) {
          double et2 = j2.Et();

          double etWeight = et1*et2/(sumEt*sumEt);
          double dPhi = deltaPhi(j1, j2);
          double cosPhi = cos(dPhi);
          if (cos(dPhi) == 1.0)  cosPhi = 0.9999;

          if (sumPt12 >  800*GeV && sumPt12 <=  850*GeV)  _hist_EEC1->fill(cosPhi, etWeight);
          if (sumPt12 >  850*GeV && sumPt12 <=  900*GeV)  _hist_EEC2->fill(cosPhi, etWeight);
          if (sumPt12 >  900*GeV && sumPt12 <= 1000*GeV)  _hist_EEC3->fill(cosPhi, etWeight);
          if (sumPt12 > 1000*GeV && sumPt12 <= 1100*GeV)  _hist_EEC4->fill(cosPhi, etWeight);
          if (sumPt12 > 1100*GeV && sumPt12 <= 1400*GeV)  _hist_EEC5->fill(cosPhi, etWeight);
          if (sumPt12 > 1400*GeV)  _hist_EEC6->fill(cosPhi, etWeight);
        }
      }
    }


    void finalize() {

      normalize(_hist_EEC1);
      normalize(_hist_EEC2);
      normalize(_hist_EEC3);
      normalize(_hist_EEC4);
      normalize(_hist_EEC5);
      normalize(_hist_EEC6);

      vector<Point2D> points1, points2, points3, points4, points5, points6;
      size_t nBins = _hist_EEC1->numBins();
      for (size_t k = 0; k < nBins/2; ++k) {

        double y1 = _hist_EEC1->bin(k+1).sumW() - _hist_EEC1->bin(nBins-(k+2)).sumW();
        double ey1 = sqrt( sqr(_hist_EEC1->bin(k+1).errW()) + sqr(_hist_EEC1->bin(nBins-(k+2)).errW()) );
        _hist_AEEC1->bin(k+1).set(y1,ey1);

        double y2 = _hist_EEC2->bin(k+1).sumW() - _hist_EEC2->bin(nBins-(k+2)).sumW();
        double ey2 = sqrt( sqr(_hist_EEC2->bin(k+1).errW()) + sqr(_hist_EEC2->bin(nBins-(k+2)).errW()) );
        _hist_AEEC2->bin(k+1).set(y2,ey2);

        double y3 = _hist_EEC3->bin(k+1).sumW() - _hist_EEC3->bin(nBins-(k+2)).sumW();
        double ey3 = sqrt( sqr(_hist_EEC3->bin(k+1).errW()) + sqr(_hist_EEC3->bin(nBins-(k+2)).errW()) );
        _hist_AEEC3->bin(k+1).set(y3,ey3);

        double y4 = _hist_EEC4->bin(k+1).sumW() - _hist_EEC4->bin(nBins-(k+2)).sumW();
        double ey4 = sqrt( sqr(_hist_EEC4->bin(k+1).errW()) + sqr(_hist_EEC4->bin(nBins-(k+2)).errW()) );
        _hist_AEEC4->bin(k+1).set(y4,ey4);

        double y5 = _hist_EEC5->bin(k+1).sumW() - _hist_EEC5->bin(nBins-(k+2)).sumW();
        double ey5 = sqrt( sqr(_hist_EEC5->bin(k+1).errW()) + sqr(_hist_EEC5->bin(nBins-(k+2)).errW()) );
        _hist_AEEC5->bin(k+1).set(y5,ey5);

        double y6 = _hist_EEC6->bin(k+1).sumW() - _hist_EEC6->bin(nBins-(k+2)).sumW();
        double ey6 = sqrt( sqr(_hist_EEC6->bin(k+1).errW()) + sqr(_hist_EEC6->bin(nBins-(k+2)).errW()) );
        _hist_AEEC6->bin(k+1).set(y6,ey6);

      }
    }


  private:

    Histo1DPtr _hist_EEC1, _hist_EEC2, _hist_EEC3, _hist_EEC4, _hist_EEC5, _hist_EEC6;
    Estimate1DPtr _hist_AEEC1, _hist_AEEC2, _hist_AEEC3, _hist_AEEC4, _hist_AEEC5, _hist_AEEC6;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2017_I1609253);


}
