// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  /// @brief Z+charm at 8 TeV
  class CMS_2017_I1634835 : public Analysis {
  public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1634835);


      /// @name Analysis methods
      ///@{

      /// Book histograms and initialise projections before the run
      void init() {
          // Initialise and register projections
          DileptonFinder zeeFinder(91.2*GeV, 0.1, Cuts::abseta < 2.1 && Cuts::pT > 20*GeV &&
                                   Cuts::abspid == PID::ELECTRON, Cuts::massIn(71*GeV, 111*GeV));
          declare(zeeFinder, "ZeeFinder");

          DileptonFinder zmumuFinder(91.2*GeV, 0.1, Cuts::abseta < 2.1 && Cuts::pT > 20*GeV &&
                                     Cuts::abspid == PID::MUON, Cuts::massIn(71*GeV, 111*GeV));
          declare(zmumuFinder, "ZmumuFinder");

          VetoedFinalState jetConstits;
          jetConstits.addVetoOnThisFinalState(zeeFinder);
          jetConstits.addVetoOnThisFinalState(zmumuFinder);

          FastJets akt04Jets(jetConstits, JetAlg::ANTIKT, 0.4);
          declare(akt04Jets, "AntiKt04Jets");

          book(_h_Z_pt_cjet, 4, 1, 1);
          book(_h_pt_cjet, 5, 1, 1);
          book(_h_Z_pt_bjet, "TMP/Z_pt_bjet", refData(4, 1, 2));
          book(_h_pt_bjet, "TMP/pt_bjet", refData(5, 1, 2));

          // book ratio histos
          book(_h_R_Z_pt_cjet, 4, 1, 2);
          book(_h_R_jet_pt_cjet, 5, 1, 2);

          book(counter_ee, "TMP/counter_ee");
          book(counter_mm, "TMP/counter_mm");
      }

      /// Perform the per-event analysis
      void analyze(const Event& event) {
          const DileptonFinder& zeeFS = apply<DileptonFinder>(event, "ZeeFinder");
          const DileptonFinder& zmumuFS = apply<DileptonFinder>(event, "ZmumuFinder");

          const Particles& zees = zeeFS.bosons();
          const Particles& zmumus = zmumuFS.bosons();

          // We did not find exactly one Z. No good.
          if (zees.size() + zmumus.size() != 1) {
              MSG_DEBUG("Did not find exactly one good Z candidate");
              vetoEvent;
          }

          Particles leptons;
          Particle zcand;
          if (zees.size() == 1) {
              leptons = zeeFS.leptons();
              zcand = zees[0];
              counter_ee->fill();
          }
          if (zmumus.size() == 1) {
              leptons = zmumuFS.leptons();
              zcand = zmumus[0];
              counter_mm->fill();
          }

          // Cluster jets
          // NB. Veto has already been applied on leptons and photons used for dressing
          const FastJets& fj = apply<FastJets>(event, "AntiKt04Jets");
          Jets goodjets = fj.jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);
          idiscardIfAnyDeltaRLess(goodjets, leptons, 0.5);

          // We don't care about events with no isolated jets
          if (goodjets.empty()) {
              MSG_DEBUG("No jets in event");
              vetoEvent;
          }

          Jets jc_final;
          Jets jb_final;

          // identification of bjets
          int n_btag = count(goodjets, hasBTag());

          for (const Jet& j : goodjets) {
              if (j.cTagged() && n_btag == 0) {
                  jc_final.push_back(j);
              }
              if (j.bTagged()) {
                  jb_final.push_back(j);
              }
          }
          // histogram filling

          if (goodjets.size() > 0) {
              if (jc_final.size() > 0) {
                  _h_pt_cjet->fill(jc_final[0].pt()/GeV);
                  _h_Z_pt_cjet->fill(zcand.pt()/GeV);
              }
              if (jb_final.size() > 0) {
                  _h_pt_bjet->fill(jb_final[0].pt()/GeV);
                  _h_Z_pt_bjet->fill(zcand.pt()/GeV);
              }
          }
      }

      /// Normalise histograms etc., after the run
      void finalize() {
          double norm = (sumOfWeights() != 0) ? crossSection() / picobarn / sumOfWeights() : 1.0;
          // account if we have electrons and muons in the sample
          if ((counter_ee->val() > 1.) && (counter_mm->val() > 1.)) {
              norm = norm / 2.;
          }

          divide(_h_pt_cjet, _h_pt_bjet, _h_R_jet_pt_cjet);
          divide(_h_Z_pt_cjet, _h_Z_pt_bjet, _h_R_Z_pt_cjet);

          scale(_h_Z_pt_cjet, norm);
          scale(_h_pt_cjet, norm);
      }

      ///@}

      /// @name Histograms

      Histo1DPtr _h_pt_cjet;
      Histo1DPtr _h_pt_bjet;
      Histo1DPtr _h_Z_pt_cjet, _h_Z_pt_bjet;

      Estimate1DPtr _h_R_jet_pt_cjet;
      Estimate1DPtr _h_R_Z_pt_cjet;

      CounterPtr counter_ee, counter_mm;
  };

  RIVET_DECLARE_PLUGIN(CMS_2017_I1634835);

}  // namespace Rivet
