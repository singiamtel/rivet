// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @brief Z production in PbPb and pp collisions at 2.76 TeV in the dimuon and dielectron decay channels
  class CMS_2014_I1322726 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2014_I1322726);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs(Cuts::abseta < 4.9);

      //Defination of ZFinder
      ZFinder zmumuFinder(fs, Cuts::abseta < 5 && Cuts::pT > 0. * GeV, PID::MUON, 60.0 * GeV, 120.0 * GeV);  //dimuon
      declare(zmumuFinder, "ZmumuFinder");

      ZFinder zeeFinder(
          fs, Cuts::abseta < 5 && Cuts::pT > 0. * GeV, PID::ELECTRON, 60.0 * GeV, 120.0 * GeV);  //dielectron
      declare(zeeFinder, "ZeeFinder");

      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      book(_h["ds/dydp-mu"], 1, 1, 1);
      book(_h["ds/dydp-el"], 2, 1, 1);
      book(_h["ds/dy-mu"], 3, 1, 1);
      book(_h["ds/dy-el"], 4, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ZFinder& zfindermu = apply<ZFinder>(event, "ZmumuFinder");  // muon
      const Particles& zmumus = zfindermu.bosons();

      const ZFinder& zfinderel = apply<ZFinder>(event, "ZeeFinder");  // electron
      const Particles& zelels = zfinderel.bosons();

      if ((zmumus.size() + zelels.size()) != 1)
        vetoEvent;

      if (zmumus.size() == 1) {                          //cout <<"zmumus " <<zmumus <<endl;
        if (abs(zmumus[0].momentum().rapidity()) < 2) {  //5 //2

          _h["ds/dydp-mu"]->fill(zmumus[0].momentum().pT() / GeV, 0.25);
          _h["ds/dy-mu"]->fill(zmumus[0].momentum().rapidity());
        }
      }

      if (zelels.size() == 1) {                             //cout <<"electron-find"<< endl;
        if (abs(zelels[0].momentum().rapidity()) < 1.44) {  //5 //1.44

          _h["ds/dydp-el"]->fill(zelels[0].momentum().pT() / GeV, 1 / 2.88);
          _h["ds/dy-el"]->fill(zelels[0].momentum().rapidity());
        }
      }
    }

    void finalize() {
      double norm = crossSection() / picobarn / sumW();

      scale(_h["ds/dydp-mu"], norm);
      scale(_h["ds/dydp-el"], norm);
      scale(_h["ds/dy-mu"], norm);
      scale(_h["ds/dy-el"], norm);
    }

    ///@}

    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}
  };

  RIVET_DECLARE_PLUGIN(CMS_2014_I1322726);

}  // namespace Rivet
