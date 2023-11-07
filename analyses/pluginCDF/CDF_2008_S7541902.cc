// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <algorithm>

namespace Rivet {


  /// @brief CDF jet pT and multiplicity distributions in W + jets events
  ///
  /// This CDF analysis provides jet pT distributions for 4 jet multiplicity bins
  /// as well as the jet multiplicity distribution in W + jets events.
  class CDF_2008_S7541902 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2008_S7541902);


    /// @name Analysis methods
    /// @{

    void init() {
      // Set up projections
      // Basic FS
      FinalState fs((Cuts::etaIn(-3.6, 3.6)));
      declare(fs, "FS");

      // Create a final state with any e-nu pair with invariant mass 65 -> 95 GeV and ET > 20 (W decay products)
      vector<pair<PdgId,PdgId> > vids;
      vids += make_pair(PID::ELECTRON, PID::NU_EBAR);
      vids += make_pair(PID::POSITRON, PID::NU_E);
      FinalState fs2(Cuts::abseta < 3.6 && Cuts::pT >= 20*GeV);
      InvMassFinalState invfs(fs2, vids, 65*GeV, 95*GeV);
      declare(invfs, "INVFS");

      // Make a final state without the W decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      declare(vfs, "VFS");
      declare(FastJets(vfs, JetAlg::CDFJETCLU, 0.4), "Jets");

      // Book histograms
      for (int i = 0 ; i < 4 ; ++i) {
        book(_histJetEt[i], i+1, 1, 1);
        book(_histJetMultRatio[i], 5, 1, i+1);
        book(_histJetMult[i], i+6, 1, 1);
        book(_numer[i],"/TMP/numer"+std::to_string(i));
        book(_denom[i],"/TMP/denom"+std::to_string(i));
      }

    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Get the W decay products (electron and neutrino)
      const InvMassFinalState& invMassFinalState = apply<InvMassFinalState>(event, "INVFS");
      const Particles&  wDecayProducts = invMassFinalState.particles();

      FourMomentum electronP, neutrinoP;
      bool gotElectron(false), gotNeutrino(false);
      for (const Particle& p : wDecayProducts) {
        FourMomentum p4 = p.momentum();
        if (p4.Et() > 20*GeV && p4.abseta() < 1.1 && p.abspid() == PID::ELECTRON) {
          electronP = p4;
          gotElectron = true;
        }
        else if (p4.Et() > 30*GeV && p.abspid() == PID::NU_E) {
          neutrinoP = p4;
          gotNeutrino = true;
        }
      }

      // Veto event if the electron or MET cuts fail
      if (!gotElectron || !gotNeutrino) vetoEvent;

      // Veto event if the MTR cut fails
      double mT2 = 2.0 * ( electronP.pT()*neutrinoP.pT() - electronP.px()*neutrinoP.px() - electronP.py()*neutrinoP.py() );
      if (sqrt(mT2) < 20*GeV ) vetoEvent;

      // Get the jets
      const JetFinder& jetProj = apply<FastJets>(event, "Jets");
      Jets theJets = jetProj.jets(cmpMomByEt, Cuts::Et > 20*GeV);
      size_t njetsA(0), njetsB(0);
      for (const Jet& j : theJets) {
        if (j.absrap() < 2.0) {
          // Fill differential histograms for top 4 jets with Et > 20
          if (njetsA < 4 && j.Et() > 20*GeV) {
            ++njetsA;
            _histJetEt[njetsA-1]->fill(j.Et()/GeV);
          }
          // Count number of jets with Et > 25 (for multiplicity histograms)
          if (j.Et() > 25*GeV)  ++njetsB;
        }
      }

      // Increment event counter
      _denom[0]->fill();

      // Jet multiplicity
      for (size_t i = 1; i <= njetsB; ++i) {
        _histJetMult[i-1]->fill(1960);
        _numer[i-1]->fill();
        if (i == 4) break;
        _denom[i]->fill();
      }
    }


    /// Finalize
    void finalize() {

      // Loop over the non-zero multiplicities
      for (size_t i = 0; i < 4; ++i) {
        if (!_denom[i]->val())  continue;
        const YODA::Estimate0D ratio = (*_numer[i]) / (*_denom[i]);
        _histJetMultRatio[i]->bin(1).set(ratio.val(), ratio.err());
      }

      // Normalize the non-ratio histograms
      for (size_t i = 0; i < 4; ++i) {
        scale(_histJetEt[i], crossSection()/picobarn/sumOfWeights());
        scale(_histJetMult[i], crossSection()/picobarn/sumOfWeights());
      }

    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _histJetEt[4];
    BinnedEstimatePtr<int> _histJetMultRatio[4];
    BinnedHistoPtr<int> _histJetMult[4];
    CounterPtr _numer[4], _denom[4];
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2008_S7541902, CDF_2008_I768579);

}
