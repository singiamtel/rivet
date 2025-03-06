// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Jet multiplicity in lepton+jets ttbar at 8 TeV
  class CMS_2016_PAS_TOP_15_006 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2016_PAS_TOP_15_006);


    /// @name Analysis methods
    /// @{

    /// Set up projections and book histograms
    void init() {

      // Complete final state
      FinalState fs;
      Cut superLooseLeptonCuts = Cuts::pt > 5*GeV;
      SpecialLeptonFinder dressedleptons(fs, superLooseLeptonCuts);
      declare(dressedleptons, "LeptonFinder");

      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressedleptons);
      declare(FastJets(fsForJets, JetAlg::ANTIKT, 0.5), "Jets");

      // Booking of histograms
      book(_normedElectronMuonHisto, "normedElectronMuonHisto", 7, 3.5, 10.5);
      book(_absXSElectronMuonHisto , "absXSElectronMuonHisto", 7, 3.5, 10.5);
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      // Select ttbar -> lepton+jets
      const SpecialLeptonFinder& dressedleptons = apply<SpecialLeptonFinder>(event, "LeptonFinder");
      vector<FourMomentum> selleptons;
      for (const DressedLepton& dressedlepton : dressedleptons.dressedLeptons()) {
        // Select good leptons
        if (dressedlepton.pT() > 30*GeV && dressedlepton.abseta() < 2.4) selleptons += dressedlepton.mom();
        // Veto loose leptons
        else if (dressedlepton.pT() > 15*GeV && dressedlepton.abseta() < 2.5) vetoEvent;
      }
      if (selleptons.size() != 1) vetoEvent;
      // Identify hardest tight lepton
      const FourMomentum lepton = selleptons[0];

      // Jets
      const FastJets& jets   = apply<FastJets>(event, "Jets");
      const Jets      jets30 = jets.jetsByPt(Cuts::pT > 30*GeV);
      int nJets = 0, nBJets = 0;
      for (const Jet& jet : jets30) {
        if (jet.abseta() > 2.5) continue;
        if (deltaR(jet.momentum(), lepton) < 0.5) continue;
        nJets += 1;
        if (jet.bTagged(Cuts::pT > 5*GeV)) nBJets += 1;
      }
      // Require >= 4 resolved jets, of which two must be b-tagged
      if (nJets < 4 || nBJets < 2) vetoEvent;

      // Fill histograms
      _normedElectronMuonHisto->fill(min(nJets, 10));
      _absXSElectronMuonHisto ->fill(min(nJets, 10));
    }


    void finalize() {
      const double ttbarXS = !std::isnan(crossSectionPerEvent()) ? crossSection() : 252.89*picobarn;
      if (std::isnan(crossSectionPerEvent()))
        MSG_INFO("No valid cross-section given, using NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV): " <<
                 ttbarXS/picobarn << " pb");

      const double xsPerWeight = ttbarXS/picobarn / sumOfWeights();
      scale(_absXSElectronMuonHisto, xsPerWeight);

      normalize(_normedElectronMuonHisto);
    }

    /// @}


    /// @brief Special dressed lepton finder
    ///
    /// Find dressed leptons by clustering all leptons and photons
    class SpecialLeptonFinder : public FinalState {
    public:

      /// Constructor
      SpecialLeptonFinder(const FinalState& fs, const Cut& cut)
        : FinalState(cut) {
        setName("CMS_2016_PAS_TOP_15_006::SpecialLeptonFinder");
        IdentifiedFinalState ifs(fs);
        ifs.acceptIdPair(PID::PHOTON);
        ifs.acceptIdPair(PID::ELECTRON);
        ifs.acceptIdPair(PID::MUON);
        declare(FastJets(ifs, JetAlg::ANTIKT, 0.1), "LeptonJets");
      }

      /// Clone on the heap
      RIVET_DEFAULT_PROJ_CLONE(SpecialLeptonFinder);

      /// Import to avoid warnings about overload-hiding
      using Projection::operator =;

      /// Retrieve the dressed leptons
      const DressedLeptons& dressedLeptons() const { return _clusteredLeptons; }

      /// Compare projections
      CmpState compare(const Projection& p) const {
        const PCmp fscmp = mkNamedPCmp(p, "LeptonJets");
        if (fscmp != CmpState::EQ) return fscmp;
        const SpecialLeptonFinder& other = dynamic_cast<const SpecialLeptonFinder&>(p);
        const bool cutcmp = _cuts == other._cuts;
        if (!cutcmp) return CmpState::NEQ;
        return CmpState::EQ;
      }

      /// Perform the calculation
      void project(const Event& e) {
        _theParticles.clear();
        _clusteredLeptons.clear();

        DressedLeptons allClusteredLeptons;
        const Jets jets = apply<FastJets>(e, "LeptonJets").jetsByPt(Cuts::pT > 5*GeV);
        for (const Jet& jet : jets) {
          Particle lepCand;
          for (const Particle& cand : jet.particles()) {
            const int absPdgId = cand.abspid();
            if (absPdgId == PID::ELECTRON || absPdgId == PID::MUON) {
              if (cand.pt() > lepCand.pt()) lepCand = cand;
            }
          }
          // Central lepton must be the major component
          if ((lepCand.pt() < jet.pt()/2.) || (!lepCand.isChargedLepton())) continue;

          DressedLepton lepton = DressedLepton(lepCand);
          for (const Particle& cand : jet.particles()) {
            if (isSame(cand, lepCand)) continue;
            lepton.addConstituent(cand, true);
          }
          allClusteredLeptons.push_back(lepton);
        }

        for (const DressedLepton& lepton : allClusteredLeptons) {
          if (_cuts->accept(static_cast<const Particle&>(lepton))) {
            _clusteredLeptons.push_back(lepton);
            _theParticles.push_back(lepton.bareLepton());
            _theParticles += lepton.photons();
          }
        }
      }

    protected:

      /// Container which stores the clustered lepton objects
      DressedLeptons _clusteredLeptons;

    };


  private:

    /// Histograms
    Histo1DPtr _normedElectronMuonHisto, _absXSElectronMuonHisto;

  };



  RIVET_DECLARE_PLUGIN(CMS_2016_PAS_TOP_15_006);

}
