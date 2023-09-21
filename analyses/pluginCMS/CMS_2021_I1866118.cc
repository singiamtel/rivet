// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Study of Z boson plus jets events using variables sensitive to double-parton scattering in pp collisions at 13 TeV
  class CMS_2021_I1866118 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1866118);

    /// Initialization
    void init() {
      const FinalState fs;

      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 27 * GeV;
      ZFinder zmumufinder(fs, cut, PID::MUON, 70 * GeV, 110 * GeV);
      declare(zmumufinder, "zmumufinder");

      // Define veto FS in order to prevent Z-decay products entering the jet algorithm

      VetoedFinalState had_fs;
      had_fs.addVetoOnThisFinalState(zmumufinder);
      FastJets jets(had_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      book(h_dphi_Z1J_cn, 1, 1, 1);
      book(h_reldpt_Z1J_cn, 2, 1, 1);
      book(h_dphi_Zdijet_Z2J_cn, 3, 1, 1);
      book(h_reldpt_Zdijet_Z2J_cn, 4, 1, 1);
      book(h_reldpt_j1j2_Z2J_cn, 5, 1, 1);

      book(_h["dphi_Z1J"], 6, 1, 1);
      book(_h["reldpt_Z1J"], 7, 1, 1);
      book(_h["dphi_Zdijet_Z2J"], 8, 1, 1);
      book(_h["reldpt_Zdijet_Z2J"], 9, 1, 1);
      book(_h["reldpt_j1j2_Z2J"], 10, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ZFinder& zmumufinder = apply<ZFinder>(event, "zmumufinder");
      const Particles& zmumus = zmumufinder.bosons();
      if (zmumus.size() != 1) {
        vetoEvent;
      }

      // Find the (dressed!) leptons
      const Particles& leptons = zmumufinder.constituents();
      if (leptons.size() != 2)
        vetoEvent;

      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20 * GeV && Cuts::abseta < 2.4);
      idiscardIfAnyDeltaRLess(jets, leptons, 0.4);

      const size_t Njets = jets.size();

      if (Njets < 1)
        vetoEvent;

      if (Njets >= 1) {
        h_dphi_Z1J_cn->fill(deltaPhi(zmumus[0], jets[0]));
        h_reldpt_Z1J_cn->fill((zmumus[0] + jets[0].momentum()).pt() / (zmumus[0].pt() + jets[0].pt()));

        _h["dphi_Z1J"]->fill(deltaPhi(zmumus[0], jets[0]));
        _h["reldpt_Z1J"]->fill((zmumus[0] + jets[0].momentum()).pt() / (zmumus[0].pt() + jets[0].pt()));
      }

      if (Njets >= 2) {
        FourMomentum dij = jets[0].momentum() + jets[1].momentum();

        h_dphi_Zdijet_Z2J_cn->fill(deltaPhi(zmumus[0], dij));
        h_reldpt_Zdijet_Z2J_cn->fill((zmumus[0] + dij).pt() / (zmumus[0].pt() + dij.pt()));
        h_reldpt_j1j2_Z2J_cn->fill(dij.pt() / (jets[0].pt() + jets[1].pt()));

        _h["dphi_Zdijet_Z2J"]->fill(deltaPhi(zmumus[0], dij));
        _h["reldpt_Zdijet_Z2J"]->fill((zmumus[0] + dij).pt() / (zmumus[0].pt() + dij.pt()));
        _h["reldpt_j1j2_Z2J"]->fill(dij.pt() / (jets[0].pt() + jets[1].pt()));
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double norm = (sumOfWeights() != 0) ? crossSection()/picobarn/sumOfWeights() : 1.0;

      scale(h_dphi_Z1J_cn, norm);
      scale(h_reldpt_Z1J_cn, norm);
      scale(h_dphi_Zdijet_Z2J_cn, norm);
      scale(h_reldpt_Zdijet_Z2J_cn, norm);
      scale(h_reldpt_j1j2_Z2J_cn, norm);

      for (auto& item : _h) {
        double rho = item.second->density(false);
        if (rho)  scale(item.second, 1.0/rho);
      }
    }

  private:
    /// @name Histogram objects
    //@{
    Histo1DPtr h_dphi_Z1J_cn;
    Histo1DPtr h_reldpt_Z1J_cn;
    Histo1DPtr h_dphi_Zdijet_Z2J_cn;
    Histo1DPtr h_reldpt_Zdijet_Z2J_cn;
    Histo1DPtr h_reldpt_j1j2_Z2J_cn;

    map<string,Histo1DPtr> _h;

    //@}
  };

  // Hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2021_I1866118);

}  // namespace Rivet
