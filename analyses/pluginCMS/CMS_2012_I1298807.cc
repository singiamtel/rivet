// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"

namespace Rivet {


  /// Inclusive ZZ production cross section and constraints on anomalous triple gauge couplings
  class CMS_2012_I1298807 : public Analysis {
  public:

    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I1298807);


    /// Initialise projections and histograms
    void init() {

      FinalState leptons((Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5) ||
                         (Cuts::abspid == PID::MUON && Cuts::abseta < 2.4));
      declare(leptons, "Leptons");

      Cut cut_el = Cuts::abseta < 2.5 && Cuts::pT > 7.0*GeV;
      DileptonFinder zeefinder(91.2*GeV, 0.1, cut_el && Cuts::abspid == PID::ELECTRON,
                               Cuts::massIn(60*GeV, 120*GeV));
      declare(zeefinder, "ZeeFinder");

      Cut cut_mu = Cuts::abseta < 2.4 && Cuts::pT > 5.0*GeV;
      DileptonFinder zmmfinder(91.2*GeV, 0.1, cut_mu && Cuts::abspid == PID::MUON,
                               Cuts::massIn(60*GeV, 120*GeV));
      declare(zmmfinder, "ZmmFinder");

      VetoedFinalState fs_woZmm;
      fs_woZmm.addVetoOnThisFinalState(zmmfinder);
      VetoedFinalState fs_woZee;
      fs_woZee.addVetoOnThisFinalState(zeefinder);

      DileptonFinder zeefinder_woZee(fs_woZee, 91.2*GeV, 0.1,
                                     cut_el && Cuts::abspid == PID::ELECTRON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zeefinder_woZee, "Zeefinder_WoZee");
      DileptonFinder zmmfinder_woZmm(fs_woZmm, 91.2*GeV, 0.1,
                                     cut_mu && Cuts::abspid == PID::MUON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zmmfinder_woZmm, "Zmmfinder_WoZmm");

      // Book histograms
      book(_h["pt_l1"],   1, 1, 1);
      book(_h["pt_z1"],   1, 1, 2);
      book(_h["pt_zz"],   1, 1, 3);
      book(_h["m_zz"],    1, 1, 4);
      book(_h["dphi_zz"], 1, 1, 5);
      book(_h["dR_zz"],   1, 1, 6);

    }


    // Perform the per-event analysis
    void analyze(const Event& evt) {

      // Find leading leptons and apply cuts
      const Particles& leptons = apply<FinalState>(evt, "Leptons").particlesByPt();
      if (leptons.size() < 2) vetoEvent;
      const double leading_l_pt = leptons[0].pT();
      const double second_l_pt = leptons[1].pT();
      if (leading_l_pt < 20*GeV || second_l_pt < 10*GeV) vetoEvent;

      // Find acceptable ZZ combinations and build four-momenta, otherwise veto
      const DileptonFinder& zeefinder = apply<DileptonFinder>(evt, "ZeeFinder");
      const DileptonFinder& zeefinder_woZee = apply<DileptonFinder>(evt, "Zeefinder_WoZee");
      const DileptonFinder& zmmfinder = apply<DileptonFinder>(evt, "ZmmFinder");
      const DileptonFinder& zmmfinder_woZmm = apply<DileptonFinder>(evt, "Zmmfinder_WoZmm");

      FourMomentum pZ_a, pZ_b, pZ_1, pZ_2;
      FourMomentum pZZ, Z_a_l1, Z_a_l2, Z_b_l1, Z_b_l2;
      if (zeefinder.bosons().size() > 0 && zmmfinder.bosons().size() > 0) {
        pZ_a = zeefinder.boson();
        pZ_b = zmmfinder.boson();
        pZZ = pZ_a + pZ_b;
        pZ_1 = pZ_a;
        pZ_2 = pZ_b;
        Z_a_l1 = zeefinder.constituents()[0];
        Z_a_l2 = zeefinder.constituents()[1];
        Z_b_l1 = zmmfinder.constituents()[0];
        Z_b_l2 = zmmfinder.constituents()[1];
      } else if (zeefinder.bosons().size() > 0 && zeefinder_woZee.bosons().size() > 0) {
        pZ_a = zeefinder.boson();
        pZ_b = zeefinder_woZee.boson();
        pZZ = pZ_a + pZ_b;
        pZ_1 = pZ_a;
        pZ_2 = pZ_b;
        Z_a_l1 = zeefinder.constituents()[0];
        Z_a_l2 = zeefinder.constituents()[1];
        Z_b_l1 = zeefinder_woZee.constituents()[0];
        Z_b_l2 = zeefinder_woZee.constituents()[1];
      } else if (zmmfinder.bosons().size() > 0 && zmmfinder_woZmm.bosons().size() > 0) {
        pZ_a = zmmfinder.boson();
        pZ_b = zmmfinder_woZmm.boson();
        pZZ = pZ_a + pZ_b;
        pZ_1 = pZ_a;
        pZ_2 = pZ_b;
        Z_a_l1 = zmmfinder.constituents()[0];
        Z_a_l2 = zmmfinder.constituents()[1];
        Z_b_l1 = zmmfinder_woZmm.constituents()[0];
        Z_b_l2 = zmmfinder_woZmm.constituents()[1];
      } else {
        vetoEvent;
      }

      // Set ordered pT variables
      /// @todo Looks like there should be a nicer way than this
      double pt_l1 = Z_a_l1.pT();
      if (Z_a_l2.pT() > pt_l1) pt_l1 = Z_a_l2.pT();
      if (Z_b_l1.pT() > pt_l1) pt_l1 = Z_b_l1.pT();
      if (Z_b_l2.pT() > pt_l1) pt_l1 = Z_b_l2.pT();

      // Leading Z pT
      double pt_z1 = pZ_a.pT();
      if (pZ_b.pT() > pZ_a.pT()) {
        pt_z1 = pZ_b.pT();
        pZ_1 = pZ_b;
        pZ_2 = pZ_a;
      }

      // Fill histograms
      _h["pt_zz"]->fill(pZZ.pT()/GeV);
      _h["m_zz"]->fill(pZZ.mass()/GeV);
      _h["dphi_zz"]->fill(deltaPhi(pZ_a, pZ_b));
      _h["dR_zz"]->fill(deltaR(pZ_a, pZ_b, PSEUDORAPIDITY));
      _h["pt_z1"]->fill(pt_z1/GeV);
      _h["pt_l1"]->fill(pt_l1/GeV);

    }


    /// Scale histograms
    /// @todo Why is it scaling by bin width twice??
    void finalize() {

      for (auto& item : _h) {
        double area = 0.;
        for (auto& b : item.second->bins()) {
          b.scaleW(1.0/b.xWidth());
          area += b.sumW() / b.dVol();
        }
        if (area)  scale(item.second, 1.0 / area);
      }
    }


    /// Histograms
    map<string, Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(CMS_2012_I1298807);

}
