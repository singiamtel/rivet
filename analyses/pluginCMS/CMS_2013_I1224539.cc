// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

namespace Rivet {


  /// CMS jet mass measurement in W, Z and dijet events
  class CMS_2013_I1224539 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2013_I1224539);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options
      WJET = true;
      ZJET = true;
      DIJET = true;
      if ( getOption("JMODE") == "W" ) {
        ZJET = false;
        DIJET = false;
      }
      if ( getOption("JMODE") == "Z" ) {
        WJET = false;
        DIJET = false;
      }
      if ( getOption("JMODE") == "DIJET" ) {
        WJET = false;
        ZJET = false;
      }

      if (WJET) {
        // filling W+jet histos

        // Find W's with pT > 120, MET > 50
        LeptonFinder lf(Cuts::abseta < 2.4 && Cuts::pT > 80*GeV && Cuts::abspid == PID::ELECTRON, 0.2);
        declare(lf, "Leptons");
        declare(MissingMomentum(), "MET");

        // W+jet jet collections
        VetoedFinalState rfs;
        rfs.vetoFinalState(lf);
        declare(FastJets(rfs, JetAlg::ANTIKT, 0.7), "JetsAK7_wj");
        declare(FastJets(rfs, JetAlg::CAM, 0.8), "JetsCA8_wj");
        declare(FastJets(rfs, JetAlg::CAM, 1.2), "JetsCA12_wj");

        // Histograms
        /// @note These are 2D histos rendered into slices
        const int wjetsOffset = 51;
        for (size_t i = 0; i < N_PT_BINS_vj; ++i) {
          book(_h_ungroomedJetMass_AK7_wj[i] ,wjetsOffset+i+1+0*N_PT_BINS_vj, 1, 1);
          book(_h_filteredJetMass_AK7_wj[i] ,wjetsOffset+i+1+1*N_PT_BINS_vj, 1, 1);
          book(_h_trimmedJetMass_AK7_wj[i] ,wjetsOffset+i+1+2*N_PT_BINS_vj, 1, 1);
          book(_h_prunedJetMass_AK7_wj[i] ,wjetsOffset+i+1+3*N_PT_BINS_vj, 1, 1);
          book(_h_prunedJetMass_CA8_wj[i] ,wjetsOffset+i+1+4*N_PT_BINS_vj, 1, 1);
          if (i > 0) book(_h_filteredJetMass_CA12_wj[i] ,wjetsOffset+i+5*N_PT_BINS_vj, 1, 1);
        }
      }

      if (ZJET) {
        // filling Z+jet histos

        // Find Zs with pT > 120 GeV
        DileptonFinder zfinder(91.2*GeV, 0.2, Cuts::abseta < 2.4 && Cuts::pT > 30*GeV &&
                               Cuts::abspid == PID::ELECTRON, Cuts::massIn(80*GeV, 100*GeV));
        declare(zfinder, "DileptonFinder");

        // Z+jet jet collections
        VetoedFinalState rfs;
        rfs.vetoFinalState(zfinder);
        declare(FastJets(rfs, JetAlg::ANTIKT, 0.7), "JetsAK7_zj");
        declare(FastJets(rfs, JetAlg::CAM, 0.8), "JetsCA8_zj");
        declare(FastJets(rfs, JetAlg::CAM, 1.2), "JetsCA12_zj");

        // Histograms
        /// @note These are 2D histos rendered into slices
        const int zjetsOffset = 28;
        for (size_t i = 0; i < N_PT_BINS_vj; ++i ) {
          book(_h_ungroomedJetMass_AK7_zj[i] ,zjetsOffset+i+1+0*N_PT_BINS_vj, 1, 1);
          book(_h_filteredJetMass_AK7_zj[i] ,zjetsOffset+i+1+1*N_PT_BINS_vj,1,1);
          book(_h_trimmedJetMass_AK7_zj[i] ,zjetsOffset+i+1+2*N_PT_BINS_vj,1,1);
          book(_h_prunedJetMass_AK7_zj[i] ,zjetsOffset+i+1+3*N_PT_BINS_vj,1,1);
          book(_h_prunedJetMass_CA8_zj[i] ,zjetsOffset+i+1+4*N_PT_BINS_vj,1,1);
          if (i > 0) book(_h_filteredJetMass_CA12_zj[i] ,zjetsOffset+i+5*N_PT_BINS_vj,1,1);
        }

      }

      if (DIJET){

        // Jet collections
        FinalState fs;
        declare(FastJets(fs, JetAlg::ANTIKT, 0.7), "JetsAK7");
        declare(FastJets(fs, JetAlg::CAM, 0.8), "JetsCA8");
        declare(FastJets(fs, JetAlg::CAM, 1.2), "JetsCA12");

        // Histograms
        for (size_t i = 0; i < N_PT_BINS_dj; ++i ) {
          book(_h_ungroomedAvgJetMass_dj[i] ,i+1+0*N_PT_BINS_dj, 1, 1);
          book(_h_filteredAvgJetMass_dj[i] ,i+1+1*N_PT_BINS_dj, 1, 1);
          book(_h_trimmedAvgJetMass_dj[i] ,i+1+2*N_PT_BINS_dj, 1, 1);
          book(_h_prunedAvgJetMass_dj[i] ,i+1+3*N_PT_BINS_dj, 1, 1);
        }

      }
    }

    bool isBackToBack_wj(const Particle& l, const P4& pmiss, const fastjet::PseudoJet& psjet) {
      const FourMomentum w = l + pmiss;
      /// @todo We should make FourMomentum know how to construct itself from a PseudoJet
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(w, jmom) > 2.0 && deltaR(l, jmom) > 1.0 && deltaPhi(pmiss, jmom) > 0.4);
    }

    bool isBackToBack_zj(const DileptonFinder& zf, const fastjet::PseudoJet& psjet) {
      const FourMomentum& z = zf.bosons()[0].momentum();
      const FourMomentum& l1 = zf.constituents()[0].momentum();
      const FourMomentum& l2 = zf.constituents()[1].momentum();
      /// @todo We should make FourMomentum know how to construct itself from a PseudoJet
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(z, jmom) > 2.0 && deltaR(l1, jmom) > 1.0 && deltaR(l2, jmom) > 1.0);
    }


    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin_vj(double ptJ) {
      const double ptBins_vj[N_PT_BINS_vj+1] = { 125.0, 150.0, 220.0, 300.0, 450.0 };
      for (size_t ibin = 0; ibin < N_PT_BINS_vj; ++ibin) {
        if (inRange(ptJ, ptBins_vj[ibin], ptBins_vj[ibin+1])) return ibin;
      }
      return N_PT_BINS_vj;
    }

    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin_jj(double ptJ) {
      const double ptBins_dj[N_PT_BINS_dj+1] = { 220.0, 300.0, 450.0, 500.0, 600.0, 800.0, 1000.0, 1500.0};
      for (size_t ibin = 0; ibin < N_PT_BINS_dj; ++ibin) {
        if (inRange(ptJ, ptBins_dj[ibin], ptBins_dj[ibin+1])) return ibin;
      }
      return N_PT_BINS_dj;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (WJET) {
        do { //< use do..while to allow breaks

          // Get the W
          /// @note Closest-matching mu+MET to mT == mW: not a great strategy!
          const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
          if (pmiss.pT() < 50*GeV) break;
          const Particles& es = apply<LeptonFinder>(event, "Leptons").particles();
          const int iefound = closestMatchIndex(es, pmiss, Kin::mT, 80.4*GeV, 50*GeV, 1000*GeV);

          // Require a fairly high-pT W and charged lepton (lepton cut already applied)
          if (iefound < 0) break;
          const Particle& e = es[iefound];
          const P4 pW = e.mom() + pmiss;
          if (pW.pT() < 120*GeV) break;

          // Get the pseudojets.
          const PseudoJets psjetsCA8_wj = apply<FastJets>(event, "JetsCA8_wj").pseudojetsByPt( 50.0*GeV );
          const PseudoJets psjetsCA12_wj = apply<FastJets>(event, "JetsCA12_wj").pseudojetsByPt( 50.0*GeV );

          // AK7 jets
          const PseudoJets psjetsAK7_wj = apply<FastJets>(event, "JetsAK7_wj").pseudojetsByPt( 50.0*GeV );
          if (!psjetsAK7_wj.empty()) {
            // Get the leading jet and make sure it's back-to-back with the W
            const fastjet::PseudoJet& j0 = psjetsAK7_wj[0];
            if (isBackToBack_wj(e, pmiss, j0)) {
              const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
              if (njetBin < N_PT_BINS_vj) {
                fastjet::PseudoJet filtered0 = _filter(j0);
                fastjet::PseudoJet trimmed0 = _trimmer(j0);
                fastjet::PseudoJet pruned0 = _pruner(j0);
                _h_ungroomedJetMass_AK7_wj[njetBin]->fill(j0.m()/GeV);
                _h_filteredJetMass_AK7_wj[njetBin]->fill(filtered0.m()/GeV);
                _h_trimmedJetMass_AK7_wj[njetBin]->fill(trimmed0.m()/GeV);
                _h_prunedJetMass_AK7_wj[njetBin]->fill(pruned0.m()/GeV);
              }
            }
          }

          // CA8 jets
          if (!psjetsCA8_wj.empty()) {
            // Get the leading jet and make sure it's back-to-back with the W
            const fastjet::PseudoJet& j0 = psjetsCA8_wj[0];
            if (isBackToBack_wj(e, pmiss, j0)) {
              const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
              if (njetBin < N_PT_BINS_vj) {
                fastjet::PseudoJet pruned0 = _pruner(j0);
                _h_prunedJetMass_CA8_wj[njetBin]->fill(pruned0.m()/GeV);
              }
            }
          }

          // CA12 jets
          if (!psjetsCA12_wj.empty()) {
            // Get the leading jet and make sure it's back-to-back with the W
            const fastjet::PseudoJet& j0 = psjetsCA12_wj[0];
            if (isBackToBack_wj(e, pmiss, j0)) {
              const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
              if (njetBin < N_PT_BINS_vj&&njetBin>0) {
                fastjet::PseudoJet filtered0 = _filter(j0);
                _h_filteredJetMass_CA12_wj[njetBin]->fill( filtered0.m() / GeV);
              }
            }
          }

          break;
        } while (true);
      }

      if (ZJET) {
        // Get the Z
        const DileptonFinder& zfinder = apply<DileptonFinder>(event, "DileptonFinder");
        if (zfinder.bosons().size() == 1) {
          const Particle& z = zfinder.bosons()[0];
          if (z.constituents().size() < 2) {
            MSG_WARNING("Found a Z with less than 2 constituents.");
            vetoEvent;
          }
          const Particle& l1 = z.constituents()[0];
          const Particle& l2 = z.constituents()[1];
          MSG_DEBUG(l1.pT() << " " << l2.pT());

          // Require a high-pT Z (and constituents)
          if (l1.pT() >= 30*GeV && l2.pT() >= 30*GeV && z.pT() >= 120*GeV) {

            // AK7 jets
            const PseudoJets& psjetsAK7_zj = apply<FastJets>(event, "JetsAK7_zj").pseudojetsByPt(50.0*GeV);
            if (!psjetsAK7_zj.empty()) {
              // Get the leading jet and make sure it's back-to-back with the Z
              const fastjet::PseudoJet& j0 = psjetsAK7_zj[0];
              if (isBackToBack_zj(zfinder, j0)) {
                const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
                if (njetBin < N_PT_BINS_vj) {
                  fastjet::PseudoJet filtered0 = _filter(j0);
                  fastjet::PseudoJet trimmed0 = _trimmer(j0);
                  fastjet::PseudoJet pruned0 = _pruner(j0);
                  _h_ungroomedJetMass_AK7_zj[njetBin]->fill(j0.m()/GeV);
                  _h_filteredJetMass_AK7_zj[njetBin]->fill(filtered0.m()/GeV);
                  _h_trimmedJetMass_AK7_zj[njetBin]->fill(trimmed0.m()/GeV);
                  _h_prunedJetMass_AK7_zj[njetBin]->fill(pruned0.m()/GeV);
                }
              }
            }

            // CA8 jets
            const PseudoJets& psjetsCA8_zj = apply<FastJets>(event, "JetsCA8_zj").pseudojetsByPt(50.0*GeV);
            if (!psjetsCA8_zj.empty()) {
              // Get the leading jet and make sure it's back-to-back with the Z
              const fastjet::PseudoJet& j0 = psjetsCA8_zj[0];
              if (isBackToBack_zj(zfinder, j0)) {
                const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
                if (njetBin < N_PT_BINS_vj) {
                  fastjet::PseudoJet pruned0 = _pruner(j0);
                  _h_prunedJetMass_CA8_zj[njetBin]->fill(pruned0.m()/GeV);
                }
              }
            }

            // CA12 jets
            const PseudoJets& psjetsCA12_zj = apply<FastJets>(event, "JetsCA12_zj").pseudojetsByPt(50.0*GeV);
            if (!psjetsCA12_zj.empty()) {
              // Get the leading jet and make sure it's back-to-back with the Z
              const fastjet::PseudoJet& j0 = psjetsCA12_zj[0];
              if (isBackToBack_zj(zfinder, j0)) {
                const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
                if (njetBin>0 && njetBin < N_PT_BINS_vj) {
                  fastjet::PseudoJet filtered0 = _filter(j0);
                  _h_filteredJetMass_CA12_zj[njetBin]->fill( filtered0.m() / GeV);
                }
              }
            }
          }
        }
      }

      if (DIJET) {
        // Look at events with >= 2 jets
        const PseudoJets& psjetsAK7 = apply<FastJets>(event, "JetsAK7").pseudojetsByPt( 50.0*GeV );
        if (psjetsAK7.size() >= 2) {

          // Get the leading two jets and find their average pT
          const fastjet::PseudoJet& j0 = psjetsAK7[0];
          const fastjet::PseudoJet& j1 = psjetsAK7[1];
          double ptAvg = 0.5 * (j0.pt() + j1.pt());

          // Find the appropriate mean pT bin and escape if needed
          const size_t njetBin = findPtBin_jj(ptAvg/GeV);
          if (njetBin < N_PT_BINS_dj) {

            // Now run the substructure algs...
            fastjet::PseudoJet filtered0 = _filter(j0);
            fastjet::PseudoJet filtered1 = _filter(j1);
            fastjet::PseudoJet trimmed0 = _trimmer(j0);
            fastjet::PseudoJet trimmed1 = _trimmer(j1);
            fastjet::PseudoJet pruned0 = _pruner(j0);
            fastjet::PseudoJet pruned1 = _pruner(j1);

            // ... and fill the histograms
            _h_ungroomedAvgJetMass_dj[njetBin]->fill(0.5*(j0.m() + j1.m())/GeV);
            _h_filteredAvgJetMass_dj[njetBin]->fill(0.5*(filtered0.m() + filtered1.m())/GeV);
            _h_trimmedAvgJetMass_dj[njetBin]->fill(0.5*(trimmed0.m() + trimmed1.m())/GeV);
            _h_prunedAvgJetMass_dj[njetBin]->fill(0.5*(pruned0.m() + pruned1.m())/GeV);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double normalizationVal = 1000;

      for (size_t i = 0; i < N_PT_BINS_vj; ++i) {
        if (WJET) {
          normalize(_h_ungroomedJetMass_AK7_wj[i], normalizationVal);
          normalize(_h_filteredJetMass_AK7_wj[i], normalizationVal);
          normalize(_h_trimmedJetMass_AK7_wj[i], normalizationVal);
          normalize(_h_prunedJetMass_AK7_wj[i], normalizationVal);
          normalize(_h_prunedJetMass_CA8_wj[i], normalizationVal);
          if (i > 0) normalize( _h_filteredJetMass_CA12_wj[i], normalizationVal);
        }
        if (ZJET) {
          normalize( _h_ungroomedJetMass_AK7_zj[i], normalizationVal);
          normalize( _h_filteredJetMass_AK7_zj[i], normalizationVal);
          normalize( _h_trimmedJetMass_AK7_zj[i], normalizationVal);
          normalize( _h_prunedJetMass_AK7_zj[i], normalizationVal);
          normalize( _h_prunedJetMass_CA8_zj[i], normalizationVal);
          if (i > 0) normalize( _h_filteredJetMass_CA12_zj[i], normalizationVal);
        }
      }
      if (DIJET) {
        for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
          normalize(_h_ungroomedAvgJetMass_dj[i], normalizationVal);
          normalize(_h_filteredAvgJetMass_dj[i], normalizationVal);
          normalize(_h_trimmedAvgJetMass_dj[i], normalizationVal);
          normalize(_h_prunedAvgJetMass_dj[i], normalizationVal);
        }
      }
    }

    /// @}


  private:

    bool WJET, ZJET, DIJET;

    /// @name FastJet grooming tools (configured in constructor init list)
    /// @{
    const fastjet::Filter _filter{fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))};
    const fastjet::Filter _trimmer{fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03))};
    const fastjet::Pruner _pruner{fastjet::Pruner(fastjet::cambridge_algorithm, 0.1, 0.5)};
    /// @}


    /// @name Histograms
    /// @{
    enum BINS_vj { PT_125_150_vj=0, PT_150_220_vj, PT_220_300_vj, PT_300_450_vj, N_PT_BINS_vj };
    // W+jet
    Histo1DPtr _h_ungroomedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_trimmedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_CA8_wj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_CA12_wj[N_PT_BINS_vj];
    //Z+jet
    Histo1DPtr _h_ungroomedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_trimmedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_CA8_zj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_CA12_zj[N_PT_BINS_vj];
    // DIJET
    enum BINS_dj { PT_220_300_dj=0, PT_300_450_dj, PT_450_500_dj, PT_500_600_dj,
                   PT_600_800_dj, PT_800_1000_dj, PT_1000_1500_dj, N_PT_BINS_dj };
    Histo1DPtr _h_ungroomedJet0pt, _h_ungroomedJet1pt;
    Histo1DPtr _h_ungroomedAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_filteredAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_trimmedAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_prunedAvgJetMass_dj[N_PT_BINS_dj];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CMS_2013_I1224539);

}
