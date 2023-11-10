// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @brief Differential Z cross section measurement in phi* at 8 TeV
  class CMS_2017_I1631985 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1631985);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // default to combined ee+mumu
      _mode = 2;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;
      if ( getOption("LMODE") == "EMU" ) _mode = 2;

      _twodim = false;
      if ( getOption("TWODIM") == "YES" ) _twodim = true;

      FinalState fs;
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;

      // The experimental result is unfolded to the 'born' level (pre-FSR), which is discouraged in Rivet
      // The implementation in Rivet here is using dressed leptons
      ZFinder zeeFind(fs, cut, PID::ELECTRON, 60.0*GeV, 120.0*GeV, 0.1, LeptonOrigin::PROMPT, PhotonOrigin::NODECAY, PhotonsAsConstituents::YES);
      declare(zeeFind, "ZeeFind");
      ZFinder zmmFind(fs, cut, PID::MUON    , 60.0*GeV, 120.0*GeV, 0.1, LeptonOrigin::PROMPT, PhotonOrigin::NODECAY, PhotonsAsConstituents::YES);
      declare(zmmFind, "ZmmFind");

      // Book histograms
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h_Zll_phiStar,      1, 1, 1);
      book(_h_Zll_phiStar_norm, 2, 1, 1);
      if (_twodim) {
        book(_h2D_Zll_phiStar_y,      3, 1, 1);
        book(_h2D_Zll_phiStar_y_norm, 4, 1, 1);
      }

      vector<double> edges = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
      book(_b_Zll_phiStar, edges);
      book(_b_Zll_phiStar_norm, edges);
      for (size_t i = 0; i < _b_Zll_phiStar->numBins(); ++i) {
        book(_b_Zll_phiStar->bin(i+1), 9, 1, i+1);
        book(_b_Zll_phiStar_norm->bin(i+1), 10, 1, i+1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zeeFS = apply<ZFinder>(event, "ZeeFind");
      const ZFinder& zmumuFS = apply<ZFinder>(event, "ZmmFind");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }
      if (zees.size() == 1 && _mode == 1)
        vetoEvent;
      if (zmumus.size() == 1 && _mode == 0)
        vetoEvent;


      bool ee_event=false;
      if( zees.size() == 1 ) ee_event = true;

      const Particles& theLeptons = ee_event ? zeeFS.constituents() : zmumuFS.constituents();
      const Particle& leadingLepton = theLeptons[0].pt() > theLeptons[1].pt() ? theLeptons[0] : theLeptons[1];

      // asymmetric cut
      if( leadingLepton.pt() < 30.0 ) vetoEvent;
      if( leadingLepton.abseta() > 2.1 ) vetoEvent;

      // calculate phi*
      const Particle& lminus = theLeptons[0].charge() < 0 ? theLeptons[0] : theLeptons[1];
      const Particle& lplus = theLeptons[0].charge() < 0 ? theLeptons[1] : theLeptons[0];
      const double thetaStar = acos(tanh( 0.5 * (lminus.eta() - lplus.eta()) ));
      const double dPhi = M_PI - deltaPhi(lminus, lplus);
      const double phiStar = tan(0.5 * dPhi) * sin(thetaStar);

      const Particle& zcand = ee_event ? zees[0] : zmumus[0];

      _h_Zll_phiStar->fill(phiStar);
      _h_Zll_phiStar_norm->fill(phiStar);

      double absRap = zcand.absrap();
      if (_twodim) {
        _h2D_Zll_phiStar_y->fill(phiStar, absRap);
        _h2D_Zll_phiStar_y_norm->fill(phiStar, absRap);
      }

      _b_Zll_phiStar->fill(absRap, phiStar);
      _b_Zll_phiStar_norm->fill(absRap, phiStar);

    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double norm = (sumOfWeights() != 0) ? crossSection()/picobarn/sumOfWeights() : 1.0;

      // when running for both ee and mm channel, need to average to get lepton xsec
      if (_mode == 2) norm /= 2.;

      scale(_h_Zll_phiStar, norm);
      scale(_b_Zll_phiStar, norm);
      divByGroupWidth(_b_Zll_phiStar);
      normalize(_h_Zll_phiStar_norm);

      if (_twodim) {
        scale(_h2D_Zll_phiStar_y, norm);
        normalize(_h2D_Zll_phiStar_y_norm);
      }

      // normalized using the sum of 2D
      normalizeGroup(_b_Zll_phiStar_norm);
      divByGroupWidth(_b_Zll_phiStar_norm);
    }

    /// @}

    Histo1DPtr _h_Zll_phiStar, _h_Zll_phiStar_norm;
    Histo1DGroupPtr _b_Zll_phiStar, _b_Zll_phiStar_norm;
    Histo2DPtr _h2D_Zll_phiStar_y, _h2D_Zll_phiStar_y_norm;

    size_t _mode;
    bool _twodim;
  };


  RIVET_DECLARE_PLUGIN(CMS_2017_I1631985);

}
