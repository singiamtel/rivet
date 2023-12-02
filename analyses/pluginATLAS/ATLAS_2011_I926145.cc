// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief Measurement of electron and muon differential cross section from heavy flavour production
  ///
  /// Lepton cross-sections differential in pT
  ///
  /// @author Paul Bell, Holger Schulz
  class ATLAS_2011_I926145 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I926145);


    /// Book histograms and initialise projections before the run
    void init() {

      // Electrons and muons
      Cut cuts = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52,  2.00)) && Cuts::pT > 7*GeV;
      IdentifiedFinalState elecs(cuts, {PID::ELECTRON, PID::POSITRON});
      declare(elecs, "elecs");
      IdentifiedFinalState muons(Cuts::abseta < 2 && Cuts::pT > 7*GeV, {PID::MUON, PID::ANTIMUON});
      declare(muons, "muons");
      IdentifiedFinalState muons_full(Cuts::abseta < 2.5 && Cuts::pT > 4*GeV, {PID::MUON, PID::ANTIMUON});
      declare(muons_full, "muons_full");

      Cut cut20 = Cuts::abseta < 2.0;
      Cut cut25 = Cuts::abseta < 2.5;

      DileptonFinder zfinder_e(91.2*GeV, 0.1, cut20 && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_e, "DileptonFinder_e");
      DileptonFinder zfinder_mu(91.2*GeV, 0.1, cut20 && Cuts::abspid == PID::MUON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_mu, "DileptonFinder_mu");
      DileptonFinder zfinder_mufull(91.2*GeV, 0.1, cut25 && Cuts::abspid == PID::MUON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_mufull, "DileptonFinder_mufull");

      LeptonFinder ef(0.2, cut20 && Cuts::abspid == PID::ELECTRON);
      declare(ef, "WElecs");
      LeptonFinder mf(0.2, cut20 && Cuts::abspid == PID::MUON);
      declare(mf, "WMuons");
      LeptonFinder mff(0.2, cut25 && Cuts::abspid == PID::MUON);
      declare(mff, "WMuonsFull");
      declare("MET", MissingMomentum());

      // Book histograms
      book(_histPt_elecs,      1, 1, 1);
      book(_histPt_muons,      2, 1, 1);
      book(_histPt_muons_full, 3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Veto event if no lepton is present
      const Particles& elecs      = apply<FinalState>(event, "elecs").particles();
      const Particles& muons      = apply<FinalState>(event, "muons").particles();
      const Particles& muons_full = apply<FinalState>(event, "muons_full").particles();
      if (elecs.empty() && muons.empty() && muons_full.empty()) vetoEvent;

      // Z veto
      const DileptonFinder& zfinder_e      = apply<DileptonFinder>(event, "DileptonFinder_e");
      const DileptonFinder& zfinder_mu     = apply<DileptonFinder>(event, "DileptonFinder_mu");
      const DileptonFinder& zfinder_mufull = apply<DileptonFinder>(event, "DileptonFinder_mufull");
      if (zfinder_e.bosons().size() > 0 || zfinder_mu.bosons().size() > 0 || zfinder_mufull.bosons().size() > 0) {
        MSG_DEBUG("Num elec Z-bosons found: " << zfinder_e.bosons().size());
        MSG_DEBUG("Num muon Z-bosons found: " << zfinder_mu.bosons().size());
        MSG_DEBUG("Num muon Z-bosons found (|eta|<2.5): " << zfinder_mufull.bosons().size());
        vetoEvent;
      }

      // W veto
      //       60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      const Particles& wes = apply<LeptonFinder>(event, "WElecs").particles();
      const int iefound = closestMatchIndex(wes, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      const Particles& wmus = apply<LeptonFinder>(event, "WMuons").particles();
      const int imfound = closestMatchIndex(wmus, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      const Particles& wfmus = apply<LeptonFinder>(event, "WMuonsFull").particles();
      const int ifmfound = closestMatchIndex(wfmus, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      if (pmiss.Et() > 25*GeV && (iefound >= 0 || imfound >= 0 || ifmfound >= 0)) {
        MSG_DEBUG("Num elec W-bosons found: " << int(iefound >= 0));
        MSG_DEBUG("Num muon W-bosons found: " << int(imfound >= 0));
        MSG_DEBUG("Num muon W-bosons found (|eta|<2.5): " << int(ifmfound >= 0));
        vetoEvent;
      }

      // Electron histogram
      if (elecs.size() > 0) {
        for (const Particle& ele : elecs) {
          if (ele.pT() < 26.0*GeV) _histPt_elecs->fill(ele.pT()/GeV);
        }
      }

      // Muon histogram
      if (muons.size() > 0) {
        for (const Particle& muo : muons) {
          if (muo.pT() < 26.0*GeV) _histPt_muons->fill(muo.pT()/GeV);
        }
      }

      // Muon full histogram
      if (muons_full.size() > 0) {
        for (const Particle& muo : muons_full) {
          if (muo.pT() < 100.0*GeV) _histPt_muons_full->fill(muo.pT()/GeV);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histPt_elecs,      crossSection()/nanobarn/sumOfWeights());
      scale(_histPt_muons,      crossSection()/nanobarn/sumOfWeights());
      scale(_histPt_muons_full, crossSection()/nanobarn/sumOfWeights());
      // the last one is microbarn on HD, but cannot sync as the d02 is
      // currently missing from HD:
      //scale(_histPt_muons_full, crossSection()/microbarn/sumOfWeights());
    }


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _histPt_elecs, _histPt_muons, _histPt_muons_full;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2011_I926145);

}
