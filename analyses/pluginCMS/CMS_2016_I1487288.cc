#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// WZ production cross section in pp collisions at 7 and 8 TeV
  class CMS_2016_I1487288 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2016_I1487288);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(Cuts::abseta < 4.9);

      FastJets fj(fs, JetAlg::ANTIKT, 0.5, JetMuons::ALL, JetInvisibles::DECAY);
      declare(fj, "Jets");

      DileptonFinder zeeFinder(91.2*GeV, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV &&
                               Cuts::abspid == PID::ELECTRON, Cuts::massIn(71*GeV, 111*GeV));
      declare(zeeFinder, "Zee");

      DileptonFinder zmumuFinder(91.2*GeV, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV &&
                                 Cuts::abspid == PID::MUON, Cuts::massIn(71*GeV, 111*GeV));
      declare(zmumuFinder, "Zmumu");

      // Initialise and register projections
      LeptonFinder ef(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON, 0.1);
      LeptonFinder mf(Cuts::abseta < 2.4 && Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON, 0.1);
      declare(ef, "Electrons");
      declare(mf, "Muons");
      declare(MissingMom(), "MET");

      book(_h_ZpT,  "d03-x01-y01");
      book(_h_Njet, "d04-x01-y01");
      book(_h_JpT,  "d05-x01-y01");

      MSG_WARNING("LIMITED VALIDITY - check info file for details!");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find Z -> l+ l-
      const DileptonFinder& zeeFS = apply<DileptonFinder>(event, "Zee");
      const DileptonFinder& zmumuFS = apply<DileptonFinder>(event, "Zmumu");
      const Particles zlls = zeeFS.bosons() + zmumuFS.bosons();
      if (zlls.empty()) vetoEvent;

      // Require some MET
      const MissingMomentum& mm = apply<MissingMom>(event, "MET");
      const P4& pmiss = mm.missingMom();
      if (pmiss.pT() < 30*GeV) vetoEvent;
	
      // Next find the W's
      const Particles& es = apply<LeptonFinder>(event, "Electrons").particles();
      const Particles es_mtfilt = select(es, [&](const Particle& e){ return inRange(mT(e, pmiss), 60*GeV, 100*GeV); });
      const int iefound = closestMatchIndex(es_mtfilt, pmiss, Kin::mass, 80.4*GeV);
      const Particles& ms = apply<LeptonFinder>(event, "Muons").particles();
      const Particles ms_mtfilt = select(ms, [&](const Particle& m){ return inRange(mT(m, pmiss), 60*GeV, 100*GeV); });
      const int imfound = closestMatchIndex(ms_mtfilt, pmiss, Kin::mass, 80.4*GeV);

      // Build a combined list of pseudo-W's
      Particles wls; wls.reserve(2);
      if (iefound < 0) wls.push_back(Particle(copysign(PID::WBOSON, es_mtfilt[iefound].pid()), es_mtfilt[iefound].mom()+pmiss));
      if (imfound < 0) wls.push_back(Particle(copysign(PID::WBOSON, ms_mtfilt[imfound].pid()), ms_mtfilt[imfound].mom()+pmiss));
      if (wls.empty()) vetoEvent;

      // If more than one Z candidate, use the one with Mll nearest to MZ
      /// @todo Can now use the closestMassIndex functions
      const Particles zlls_mz = sortBy(zlls, [](const Particle& a, const Particle& b){
          return fabs(a.mass() - 91.2*GeV) < fabs(b.mass() - 91.2*GeV); });
      const Particle& Z = zlls_mz.front();
      // const bool isZee = any(Z.constituents(), hasAbsPID(PID::ELECTRON));

      // If more than one W candidate, use the one with Mlv nearest to MW
      /// @todo Can now use the closestMassIndex functions
      const Particles wls_mw = sortBy(wls, [](const Particle& a, const Particle& b){
          return fabs(a.mass() - 80.4*GeV) < fabs(b.mass() - 80.4*GeV); });
      const Particle& W = wls_mw.front();
      // const bool isWe = any(W.constituents(), hasAbsPID(PID::ELECTRON));

      // Isolate W and Z charged leptons from each other
      for (const Particle& lw : W.constituents()) {
        if (lw.charge3() == 0) continue;
        for (const Particle& lz : Z.constituents()) {
          if (deltaR(lw, lz) < 0.1) vetoEvent;
        }
      }

      // Fill Z pT histogram
      _h_ZpT->fill(Z.pT()/GeV);


      // Isolate jets from W and Z charged leptons
      const Particles wzleps = select(W.constituents()+Z.constituents(), isChargedLepton);
      const Jets& jets = apply<FastJets>("Jets", event).jetsByPt(Cuts::pT > 30*GeV and Cuts::abseta < 2.5);
      const Jets isojets = discardIfAnyDeltaRLess(jets, wzleps, 0.5);

      // Fill jet histograms
      _h_Njet->fill(isojets.size());
      if (!isojets.empty()) _h_JpT->fill(isojets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Total cross-section is corrected for BR(W->lnu), BR(Z->ll), leptonic-tau fraction f_tau = 6.5-7.6%,
      // and unpublished detector/acceptance signal efficiencies epsilon_sig. Fix to published value: valid for shape comparison only
      const double xsec8tev = 24.09; // picobarn;
      normalize(_h_ZpT,  xsec8tev);
      normalize(_h_Njet, xsec8tev);
      normalize(_h_JpT,  xsec8tev);
    }


  private:

    /// Histogram
    Histo1DPtr _h_ZpT, _h_JpT;
    BinnedHistoPtr<int> _h_Njet;

  };



  RIVET_DECLARE_PLUGIN(CMS_2016_I1487288);

}
