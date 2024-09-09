// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// @brief ZZ production at 13.6 TeV
  class ATLAS_2023_I2723369 : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2723369);

    void init() {

      /// Dressed leptons
      Cut cut_lep = (Cuts::abspid == PID::MUON && Cuts::abseta < 2.5 && Cuts::pT > 5*GeV) ||
                    (Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && Cuts::pT > 7*GeV);
	    PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState leptons(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
	    LeptonFinder dLeptons(leptons, photons, 0.1, cut_lep);
      declare(dLeptons, "AllLeptons");

      /// Jet inputs
      FinalState fs_jet(Cuts::abseta < 5.0);
      VetoedFinalState jet_input(fs_jet);

      // reject all leptons dressed with only prompt photons from jet input
      FinalState all_leptons(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
   		LeptonFinder reject_leptons(all_leptons, photons, 0.1);
      jet_input.addVetoOnThisFinalState(reject_leptons);

      // reject prompt invisibles, including from tau decays
      VetoedFinalState invis_fs_jet(fs_jet);
      invis_fs_jet.addVetoOnThisFinalState(VisibleFinalState(fs_jet));
      PromptFinalState invis_pfs_jet(invis_fs_jet, TauDecaysAs::PROMPT);
      jet_input.addVetoOnThisFinalState(invis_pfs_jet);

      // declare jets
      FastJets jets(jet_input, JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::DECAY);
      declare(jets, "Jets");

      // Book histograms
      book(_h["m4l"],  1, 1, 1);
      book(_h["pt4l"], 2, 1, 1);

    }

    /// Do the per-event analysis
    void analyze(const Event& e) {


      const DressedLeptons& all_leps = apply<LeptonFinder>(e, "AllLeptons").dressedLeptons();
      size_t n_parts = all_leps.size();
      size_t n_OSSF_pairs = 0;

      // Form Z candidate (opposite-sign same-flavour) lepton pairs
      std::vector<Zstate> dileptons;
      for (size_t i = 0; i < n_parts; ++i) {
        for (size_t j = i + 1; j < n_parts; ++j) {
          if (isOSSF(all_leps[i], all_leps[j])){
            n_OSSF_pairs += 1;
            // Set positive charge first for later ME calculation
            if (all_leps[i].charge() > 0) {
              dileptons.emplace_back(make_pair(all_leps[i], all_leps[j]));
            } else {
              dileptons.emplace_back(make_pair(all_leps[j], all_leps[i]));
            }
          }
        }
      }

      // At least two pairs required to select ZZ->llll final state
      if (n_OSSF_pairs < 2) vetoEvent;

      // Form the quadruplet of two lepon pairs passing kinematics cuts
      std::vector<Quadruplet> quadruplets;
      for (size_t i = 0; i < dileptons.size(); ++i) {
        for (size_t j = i+1; j < dileptons.size(); ++j) {
          // Only use unique leptons
          if (isSame( dileptons[i].first , dileptons[j].first  )) continue;
          if (isSame( dileptons[i].first , dileptons[j].second )) continue;
          if (isSame( dileptons[i].second, dileptons[j].first  )) continue;
          if (isSame( dileptons[i].second, dileptons[j].second )) continue;

          Particles leptons{ dileptons[i].first, dileptons[i].second,
                             dileptons[j].first, dileptons[j].second };
          isortByPt(leptons);

          // Apply kinematic cuts
          if (leptons[0].pT() < 20*GeV)  continue;
          if (leptons[1].pT() < 10*GeV)  continue;

          // Form the quad with pair closest to Z pole first
          if (dileptons[i].Zdist() < dileptons[j].Zdist()) {
            quadruplets.emplace_back(dileptons[i], dileptons[j]);
          } else {
            quadruplets.emplace_back(dileptons[j], dileptons[i]);
          }
        }
      }

      // Veto if no quad passes kinematic selection
      size_t n_quads = quadruplets.size();
      if (n_quads == 0) vetoEvent;

      // To resolve ambiguities in lepton pairing order quads by channel priority first, then m12 - mz and m34 - mz
      // The first in every channel is considered nominal
      std::sort(quadruplets.begin(), quadruplets.end(),
        [](const Quadruplet& q1, const Quadruplet& q2) {
          if (q1.ch_priority == q2.ch_priority) {
            // if rarely, Z1 the same distance from the Z pole, compare Z2
            if (fabs( q1.Z1().Zdist() - q2.Z1().Zdist() ) < 1.e-5) {
              return q1.Z2().Zdist() < q2.Z2().Zdist();
            }
            return q1.Z1().Zdist() < q2.Z1().Zdist();
          }
          return q1.ch_priority < q2.ch_priority;
      });


      // Select the best quad
      Particles leptons_sel4l;
      Quadruplet quadSel;
      bool atleastonequadpassed = false;

      for (size_t iquad = 0; iquad < n_quads; ++iquad) {

        // Veto event if nominal quad was not selected in 4 lepton case
        if (n_parts == 4 && iquad > 0) vetoEvent;

        const Quadruplet& quad = quadruplets[iquad];

        // Z invariant mass requirements
        if (!(inRange(quad.Z1().mom().mass(), 66*GeV, 116*GeV))) continue;
        if (!(inRange(quad.Z2().mom().mass(), 66*GeV, 116*GeV))) continue;

        // Lepton separation and J/Psi veto
        bool b_pass_leptonseparation = true;
        bool b_pass_jpsi = true;
        leptons_sel4l.clear();
        leptons_sel4l.push_back(quad.Z1().first);
        leptons_sel4l.push_back(quad.Z1().second);
        leptons_sel4l.push_back(quad.Z2().first);
        leptons_sel4l.push_back(quad.Z2().second);

        for (size_t i = 0; i < 4; ++i) {
          for (size_t j = i+1; j < 4; ++j) {
            if ( deltaR( leptons_sel4l[i], leptons_sel4l[j]) < 0.1) b_pass_leptonseparation = false;
            if ( isOSSF(leptons_sel4l[i], leptons_sel4l[j]) && \
                 (leptons_sel4l[i].mom() + leptons_sel4l[j].mom()).mass() <= 5.*GeV) b_pass_jpsi = false;
          }
        }
        if (b_pass_leptonseparation == false || b_pass_jpsi == false) continue;

        isortByPt(leptons_sel4l);
        atleastonequadpassed = true;
        quadSel = quad;
        break;

      }

      if (!atleastonequadpassed)  vetoEvent;

      // Veto if quad not in Higgs mass window
      const FourMomentum ZZ = quadSel.mom();
      _h["pt4l"]->fill(ZZ.pT()/GeV);
      _h["m4l"]->fill(ZZ.mass()/GeV);

    }


    void finalize() {
      scale(_h, crossSection() / femtobarn  / sumOfWeights());
    }

  private:

    map<string, Histo1DPtr> _h;

    /// Generic Z candidate
    struct Zstate : public ParticlePair {

      Zstate() { }

      Zstate(ParticlePair&& _particlepair) : ParticlePair(std::move(_particlepair)) { }

      FourMomentum mom() const { return first.momentum() + second.momentum(); }

      double Zdist() const { return fabs(mom().mass() -  91.1876*GeV); }

      int flavour() const { return first.abspid(); }
    };

    /// Generic quadruplet
    struct Quadruplet {

      // find out which type it is: 4mu = 0, 4e = 1, 2mu2e = 2, 2e2mu = 3 (mm, ee, me, em)
      // channel priority is 4m, 2e2m, 2m2e, 4e
      enum class FlavCombi { mm=0, ee, me, em, undefined };

      Zstate _z1, _z2;

      FlavCombi _type;

      int ch_priority;

      Quadruplet() { }

      Quadruplet(const Zstate& z1, const Zstate& z2) : _z1(z1), _z2(z2) {
        if (     _z1.flavour() == 13 && _z2.flavour() == 13) { _type = FlavCombi::mm; ch_priority = 0; }
        else if (_z1.flavour() == 11 && _z2.flavour() == 11) { _type = FlavCombi::ee; ch_priority = 3; }
        else if (_z1.flavour() == 13 && _z2.flavour() == 11) { _type = FlavCombi::me; ch_priority = 2; }
        else if (_z1.flavour() == 11 && _z2.flavour() == 13) { _type = FlavCombi::em; ch_priority = 1; }
        else  { _type = FlavCombi::undefined; }
      }

      const Zstate& Z1() const { return _z1; }

      const Zstate& Z2() const { return _z2; }

      FourMomentum mom() const { return _z1.mom() + _z2.mom(); }

      const FlavCombi& type() const { return _type; }
    };


  };

  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2723369);
}
