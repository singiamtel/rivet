// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// ATLAS W pT analysis at 7 TeV
  class ATLAS_2011_I925932 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I925932);


    /// @name Analysis methods
    /// @{

    void init() {
      // Set up projections
      Cut cuts = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      declare("MET", MissingMomentum());
      LeptonFinder ef_dressed(0.2, cuts && Cuts::abspid == PID::ELECTRON);
      declare(ef_dressed, "Elecs_dressed");
      LeptonFinder ef_bare(0.0, cuts && Cuts::abspid == PID::ELECTRON);
      declare(ef_bare, "Elecs_bare");
      LeptonFinder mf_dressed(0.2, cuts && Cuts::abspid == PID::MUON);
      declare(mf_dressed, "Muons_dressed");
      LeptonFinder mf_bare(0.0, cuts && Cuts::abspid == PID::MUON);
      declare(mf_bare, "Muons_bare");

      // Book histograms
      book(_hist_wpt_dressed_el, 1, 1, 2);
      book(_hist_wpt_bare_el,    1, 1, 3);
      book(_hist_wpt_dressed_mu, 2, 1, 2);
      book(_hist_wpt_bare_mu,    2, 1, 3);
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // W reco, starting with MET
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.Et() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& esd = apply<LeptonFinder>(event, "Elecs_dressed").particles();
      const int iedfound = closestMatchIndex(esd, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);
      const Particles& esb = apply<LeptonFinder>(event, "Elecs_bare").particles();
      const int iebfound = closestMatchIndex(esb, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);
      const Particles& musd = apply<LeptonFinder>(event, "Muons_dressed").particles();
      const int imdfound = closestMatchIndex(musd, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);
      const Particles& musb = apply<LeptonFinder>(event, "Muons_bare").particles();
      const int imbfound = closestMatchIndex(musb, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);

      MSG_DEBUG("Found " << int(iedfound >= 0) + int(imdfound >= 0) << " dressed W -> e/mu nu");
      MSG_DEBUG("Found " << int(iebfound >= 0) + int(imbfound >= 0) << " bare W -> e/mu nu");
      if (iedfound < 0 && iebfound < 0 && imdfound < 0 && imbfound < 0) {
        MSG_DEBUG("No W bosons found");
        vetoEvent;
      }

      // "Dressed" electron
      if (iedfound >= 0) {
        const Particle& l = esd[iedfound];
        if (mT(pmiss, l) > 40*GeV) _hist_wpt_dressed_el->fill((l.mom() + pmiss).pT()/GeV);
      }

      // "Bare" electron
      if (iebfound >= 0) {
        const Particle& l = esb[iebfound];
        if (mT(pmiss, l) > 40*GeV) _hist_wpt_bare_el->fill((l.mom() + pmiss).pT()/GeV);
      }

      // "Dressed" muon
      if (imdfound >= 0) {
        const Particle& l = musd[imdfound];
        if (mT(pmiss, l) > 40*GeV) _hist_wpt_dressed_mu->fill((l.mom() + pmiss).pT()/GeV);
      }

      // "Bare" muon
      if (imbfound >= 0) {
        const Particle& l = musb[imbfound];
        if (mT(pmiss, l) > 40*GeV) _hist_wpt_bare_mu->fill((l.mom() + pmiss).pT()/GeV);
      }

    }


    // Normalize histos
    void finalize() {
      normalize(_hist_wpt_dressed_el);
      normalize(_hist_wpt_bare_el);
      normalize(_hist_wpt_dressed_mu);
      normalize(_hist_wpt_bare_mu);
    }

    /// @}


  private:

    Histo1DPtr _hist_wpt_dressed_el, _hist_wpt_bare_el, _hist_wpt_dressed_mu, _hist_wpt_bare_mu;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2011_I925932);

}
