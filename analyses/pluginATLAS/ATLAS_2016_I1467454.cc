// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief High-mass Drell-Yan at 8 TeV
  class ATLAS_2016_I1467454 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1467454);
    /// @}


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;

      const FinalState fs;
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 30*GeV;
      ZFinder zfinder(fs, cuts, _mode? PID::MUON : PID::ELECTRON, 116*GeV, 1500*GeV, 0.1);
      declare(zfinder, "ZFinder");

      size_t ch = _mode? 11 : 0; // offset
      book(_hist_mll, 18 + ch, 1, 1);

      const vector<double> mll_bins{ 116., 150., 200., 300., 500., 1500. };
      book(_hist_rap, mll_bins);
      book(_hist_deta, mll_bins);
      for (size_t i=0; i < _hist_rap->numBins(); ++i) {
        book(_hist_rap->bin(i+1),  19 + ch + i, 1, 1);
        book(_hist_deta->bin(i+1), 24 + ch + i, 1, 1);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1)  vetoEvent;

      const Particle z0  = zfinder.bosons()[0];
      /// @todo Could use z0.constituents()
      const Particle el1 = zfinder.constituentLeptons()[0];
      const Particle el2 = zfinder.constituentLeptons()[1];

      if (el1.pT() > 40*GeV || el2.pT() > 40*GeV) {
        const double mass = z0.mass();
        _hist_mll->fill(mass/GeV);
        _hist_rap->fill(mass/GeV, z0.absrap());
        _hist_deta->fill(mass/GeV, deltaEta(el1,el2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection()/sumOfWeights();
      scale(_hist_mll, sf);
      scale(_hist_rap, sf*0.5);
      scale(_hist_deta, sf*0.5);

    }

    /// @}


    /// Choose to work in electron or muon mode
    size_t _mode;


    /// @name Histograms
    /// @{
    Histo1DPtr _hist_mll;
    Histo1DGroupPtr _hist_rap, _hist_deta;
    /// @}

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1467454);

}
