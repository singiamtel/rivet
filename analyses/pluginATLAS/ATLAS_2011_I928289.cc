// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// W inclusive cross-sections at 7 TeV
  class ATLAS_2011_I928289 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I928289);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 0;
      if ( getOption("MODE") == "Z" ) _mode = 1;
      else if ( getOption("MODE") == "W" ) _mode = 2;

      // Initialise and register projections

      Cut cut = Cuts::pT >= 20.0*GeV;
      DileptonFinder zfinder_ee_bare(   91.2*GeV, 0.0, cut && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
      DileptonFinder zfinder_ee_dressed(91.2*GeV, 0.1, cut && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
      DileptonFinder zfinder_mm_bare(   91.2*GeV, 0.0, cut && Cuts::abspid == PID::MUON    , Cuts::massIn(66*GeV, 116*GeV));
      DileptonFinder zfinder_mm_dressed(91.2*GeV, 0.1, cut && Cuts::abspid == PID::MUON    , Cuts::massIn(66*GeV, 116*GeV));

      declare(zfinder_ee_bare   , "ee_bare"   );
      declare(zfinder_ee_dressed, "ee_dressed");
      declare(zfinder_mm_bare   , "mm_bare"   );
      declare(zfinder_mm_dressed, "mm_dressed");

      declare("MET", MissingMomentum());
      LeptonFinder ef_bare(0.0, Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON);
      LeptonFinder ef_dressed(0.1, Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON);
      LeptonFinder mf_bare(0.0, Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON);
      LeptonFinder mf_dressed(0.1, Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON);
      declare(ef_bare, "elecs_bare");
      declare(ef_dressed, "elecs_dressed");
      declare(mf_bare, "muons_bare");
      declare(mf_dressed, "muons_dressed");

      /// Book histograms here
      if (_mode == 0 || _mode == 1) {
        book(_h["ee_bare"],    1, 1, 1);
        book(_h["ee_dressed"], 1, 1, 2);
        book(_h["mm_bare"],    1, 1, 3);
        book(_h["mm_dressed"], 1, 1, 4);
      }
      if (_mode == 0 || _mode == 2) {
        book(_h["minus_el_bare"],    3, 1, 1);
        book(_h["minus_el_dressed"], 3, 1, 2);
        book(_h["minus_mu_bare"],    3, 1, 3);
        book(_h["minus_mu_dressed"], 3, 1, 4);
        book(_h["plus_el_bare"],     5, 1, 1);
        book(_h["plus_el_dressed"],  5, 1, 2);
        book(_h["plus_mu_bare"],     5, 1, 3);
        book(_h["plus_mu_dressed"],  5, 1, 4);
        book(_e["el_bare"],          7, 1, 1);
        book(_e["el_dressed"],       7, 1, 2);
        book(_e["mu_bare"],          7, 1, 3);
        book(_e["mu_dressed"],       7, 1, 4);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_mode == 0 || _mode == 1) {
        fillHistosZ(event);
      }
      if (_mode == 0 || _mode == 2) {
        fillHistosW(event);
      }

    }

    void fillHistosW(const Event& event) {
      // W reco, starting with MET
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.Et() < 25*GeV)  return;

      const Particles& esb = apply<LeptonFinder>(event, "elecs_bare").particles();
      const Particles esb_mtfilt = select(esb, [&](const Particle& e){ return mT(e, pmiss) > 40*GeV; });
      const int iebfound = closestMatchIndex(esb_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (iebfound >= 0) {
        const Particle& l = esb_mtfilt[iebfound];
        (l.charge3() > 0 ? _h["plus_el_bare"] : _h["minus_el_bare"])->fill(l.abseta());
      }

      const Particles& esd = apply<LeptonFinder>(event, "elecs_dressed").particles();
      const Particles esd_mtfilt = select(esd, [&](const Particle& e){ return mT(e, pmiss) > 40*GeV; });
      const int iedfound = closestMatchIndex(esd_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (iedfound >= 0) {
        const Particle& l = esd_mtfilt[iedfound];
        (l.charge3() > 0 ? _h["plus_el_dressed"] : _h["minus_el_dressed"])->fill(l.abseta());
      }

      const Particles& musb = apply<LeptonFinder>(event, "muons_bare").particles();
      const Particles musb_mtfilt = select(musb, [&](const Particle& m){ return mT(m, pmiss) > 40*GeV; });
      const int imbfound = closestMatchIndex(musb_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (imbfound >= 0) {
        const Particle& l = musb_mtfilt[imbfound];
        (l.charge3() > 0 ? _h["plus_mu_bare"] : _h["minus_mu_bare"])->fill(l.abseta());
      }

      const Particles& musd = apply<LeptonFinder>(event, "muons_dressed").particles();
      const Particles musd_mtfilt = select(musd, [&](const Particle& m){ return mT(m, pmiss) > 40*GeV; });
      const int imdfound = closestMatchIndex(musd_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (imdfound >= 0) {
        const Particle& l = musd_mtfilt[imdfound];
        (l.charge3() > 0 ? _h["plus_mu_dressed"] : _h["minus_mu_dressed"])->fill(l.abseta());
      }
    }


    void fillHistosZ(const Event& event) {

      const DileptonFinder& zfinder_ee_bare     = apply<DileptonFinder>(event, "ee_bare"   );
      const DileptonFinder& zfinder_ee_dressed  = apply<DileptonFinder>(event, "ee_dressed");
      const DileptonFinder& zfinder_mm_bare     = apply<DileptonFinder>(event, "mm_bare"   );
      const DileptonFinder& zfinder_mm_dressed  = apply<DileptonFinder>(event, "mm_dressed");

      fillPlots1D(zfinder_ee_bare,    "ee_bare");
      fillPlots1D(zfinder_ee_dressed, "ee_dressed");
      fillPlots1D(zfinder_mm_bare,    "mm_bare");
      fillPlots1D(zfinder_mm_dressed, "mm_dressed");

    }


    void fillPlots1D(const DileptonFinder& zfinder, const string& label) {
      if (zfinder.bosons().size() != 1) return;
      const FourMomentum zmom = zfinder.bosons()[0].momentum();
      _h[label]->fill(zmom.absrap());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Print summary info
      const double xs_pb = crossSection() / picobarn;
      const double sumw = sumOfWeights();
      const double sf = 0.5 * xs_pb / sumw; // 0.5 accounts for rapidity bin width

      scale(_h, sf);
      if (_mode == 0 || _mode == 2) {
        // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta)
        for (const string& type : vector<string>{"el_bare", "el_dressed", "mu_bare", "mu_dressed"}) {
          asymm(_h["plus_"+type], _h["minus_"+type], _e[type]);
        }
      }
    }

    /// @}


  private:

    size_t _mode;

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Estimate1DPtr> _e;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2011_I928289);

}
