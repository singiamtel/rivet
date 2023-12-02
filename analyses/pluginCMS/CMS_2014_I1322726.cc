// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// Z production in PbPb and pp collisions at 2.76 TeV in the dimuon and dielectron decay channels
  class CMS_2014_I1322726 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2014_I1322726);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      DileptonFinder zmumuFinder(91.2*GeV, 0.1, Cuts::abseta < 4.9 && Cuts::abspid == PID::MUON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zmumuFinder, "ZmumuFinder");

      DileptonFinder zeeFinder(91.2*GeV, 0.1, Cuts::abseta < 4.9 && Cuts::abspid == PID::ELECTRON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zeeFinder, "ZeeFinder");

      FastJets jetfs(FinalState(Cuts::abseta < 4.9), JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::NONE);
      declare(jetfs, "jets");

      book(_h["ds/dydp-mu"], 1, 1, 1);
      book(_h["ds/dydp-el"], 2, 1, 1);
      book(_h["ds/dy-mu"], 3, 1, 1);
      book(_h["ds/dy-el"], 4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const DileptonFinder& zfindermu = apply<DileptonFinder>(event, "ZmumuFinder");  // muon
      const Particles& zmumus = zfindermu.bosons();

      const DileptonFinder& zfinderel = apply<DileptonFinder>(event, "ZeeFinder");  // electron
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


    /// Histograms
    map<string, Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(CMS_2014_I1322726);

}
