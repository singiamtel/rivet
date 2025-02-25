// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/DISDiffHadron.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief H1 diffractive dijets
  ///
  /// Diffractive dijets H1 with 920 GeV p and 27.5 GeV e
  /// Tagged protons & jets found in gamma*p rest frame.
  ///
  /// @author Christine O. Rasmussen
  class H1_2015_I1343110 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2015_I1343110);

    /// @name Analysis methods
    /// @{

    // Book projections and histograms
    void init() {

      declare(DISKinematics(), "Kinematics");
      declare(DISFinalState(DISFrame::LAB), "FSLab");
      declare(FastJets(DISFinalState(DISFrame::HCM),
        fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 1.0,
        JetMuons::ALL, JetInvisibles::NONE, nullptr), "DISFSJets");
      declare(DISDiffHadron(), "Hadron");

      // Book histograms from REF data
      book(_h_PHO_sig_sqrts, 1, 1, 1);
      book(_h_DIS_sig_sqrts, 2, 1, 1);
      book(_h_PHODIS_sqrts, 3, 1, 1);

      book(_h_DIS_dsigdz, 4, 1, 1);
      book(_h_DIS_dsigdxPom, 5, 1, 1);
      book(_h_DIS_dsigdy, 6, 1, 1);
      book(_h_DIS_dsigdQ2, 7, 1, 1);
      book(_h_DIS_dsigdEtj1, 8, 1, 1);
      book(_h_DIS_dsigdMX, 9, 1, 1);
      book(_h_DIS_dsigdDeltaEta, 10, 1, 1);
      book(_h_DIS_dsigdAvgEta, 11, 1, 1);

      book(_h_PHO_dsigdz, 12, 1, 1);
      book(_h_PHO_dsigdxPom, 13, 1, 1);
      book(_h_PHO_dsigdy, 14, 1, 1);
      book(_h_PHO_dsigdxGam, 15, 1, 1);
      book(_h_PHO_dsigdEtj1, 16, 1, 1);
      book(_h_PHO_dsigdMX, 17, 1, 1);
      book(_h_PHO_dsigdDeltaEta, 18, 1, 1);
      book(_h_PHO_dsigdAvgEta, 19, 1, 1);

      book(_h_PHODIS_deltaEta, 20, 1, 1);
      book(_h_PHODIS_y, 21, 1, 1);
      book(_h_PHODIS_z, 22, 1, 1);
      book(_h_PHODIS_Etj1, 23, 1, 1);

      isPHO  = false;
      nVeto0 = 0;
      nVeto1 = 0;
      nVeto2 = 0;
      nVeto3 = 0;
      nVeto4 = 0;
      nVeto5 = 0;
      nVeto6 = 0;
      nPHO   = 0;
      nDIS   = 0;
    }

    // Do the analysis
    void analyze(const Event& event) {

      // Event weight
      isPHO  = false;

      // Projections - special handling of events where no proton found:
      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");
      if (kin.failed()) vetoEvent;
      Particle hadronOut;
      Particle hadronIn;
      try {
        const DISDiffHadron& diffhadr = apply<DISDiffHadron>(event, "Hadron");
        hadronOut = diffhadr.out();
        hadronIn  = diffhadr.in();
      } catch (const Error& e){
        vetoEvent;
      }
      ++nVeto0;

      // Determine kinematics: H1 has +z = proton direction
      int dir   = kin.orientation();
      double y  = kin.y();
      double Q2 = kin.Q2();

      // Separate into DIS and PHO regimes else veto
      if (Q2 < 2.*GeV2 && inRange(y, 0.2, 0.70)) {
        isPHO = true;
        ++nPHO;
      } else if (inRange(Q2, 4.0*GeV2, 80.*GeV2) && inRange(y, 0.2, 0.7)) {
        isPHO = false;
        ++nDIS;
      } else vetoEvent;
      ++nVeto1;

      // Find diffractive variables as defined in paper.
      // Note tagged protons in VFPS => smaller allowed xPom range
      // xPom = 1 - E'/E, t = (P-P')^2
      const double abst = abs((hadronIn.mom() - hadronOut.mom()).mass2());
      const double xPom = 1. - hadronOut.energy() / hadronIn.energy();
      // Veto if outside allowed region
      if (abst > 0.6 * GeV2) vetoEvent;
      ++nVeto2;
      if (!inRange(xPom, 0.010, 0.024)) vetoEvent;
      ++nVeto3;

      // Jet selection. Note jets are found in HCM frame, but
      // eta cut is applied in lab frame!
      Jets jets   = apply<FastJets>(event, "DISFSJets").jets(Cuts::Et > 4.* GeV, cmpMomByEt);
      // Veto if not dijets and if Et_j1 < 5.5
      if (jets.size() < 2)          vetoEvent;
      if (jets[0].Et() < 5.5 * GeV) vetoEvent;
      ++nVeto4;

      const double EtJet1 = jets[0].Et();

      // Transform from HCM to LAB and apply eta cut
      const LorentzTransform boost = kin.boostHCM();
      for (int i = 0; i < 2; ++i) {
        jets[i].transformBy(boost.inverse());
        if (!inRange(dir * jets[i].eta(),-1,2.5)) vetoEvent;
      }
      ++nVeto5;

      double sumE(0.), sumPz(0.);
      FourMomentum tot_mom;
      for(const Particle & part : apply<FinalState>(event,"FSLab").particles()) {
        if (part.genParticle() != hadronOut.genParticle()) {
          sumE += part.energy();
          sumPz += part.pz();
          tot_mom += part.mom();
        }
      }
      double EplusPz = sumE + dir * sumPz;
      double EminusPz = sumE - dir * sumPz;
      double M2X = tot_mom.mass2();

      // Pseudorapidity distributions are examined in lab frame:
      double deltaEtaJets = abs(dir * jets[0].eta() - dir * jets[1].eta());
      double avgEtaJets   = 0.5 * (dir * jets[0].eta() + dir * jets[1].eta());

      // Evaluate observables
      double zPomJets(0.), xGamJets(0.);
      FourMomentum pj = jets[0].momentum()+jets[1].momentum();
      if (isPHO){
        zPomJets = (pj.E() + dir * pj.pz()) / EplusPz;
        xGamJets = (pj.E() - dir * pj.pz()) / EminusPz;
      } else {
        zPomJets = (Q2 + pj.mass2()) / (Q2 + M2X);
      }
      // Veto events with zPom > 0.8
      if (zPomJets > 0.8) vetoEvent;
      ++nVeto6;

      // Now fill histograms
      if (isPHO){
        _h_PHO_sig_sqrts     ->fill(sqrtS()/GeV);
        _h_PHO_dsigdz        ->fill(zPomJets);
        _h_PHO_dsigdxPom     ->fill(xPom);
        _h_PHO_dsigdy        ->fill(y);
        _h_PHO_dsigdxGam     ->fill(xGamJets);
        _h_PHO_dsigdEtj1     ->fill(EtJet1);
        _h_PHO_dsigdMX       ->fill(sqrt(M2X)/GeV);
        _h_PHO_dsigdDeltaEta ->fill(deltaEtaJets);
        _h_PHO_dsigdAvgEta   ->fill(avgEtaJets);
      } else {
        _h_DIS_sig_sqrts     ->fill(sqrtS()/GeV);
        _h_DIS_dsigdz        ->fill(zPomJets);
        _h_DIS_dsigdxPom     ->fill(xPom);
        _h_DIS_dsigdy        ->fill(y);
        _h_DIS_dsigdQ2       ->fill(Q2);
        _h_DIS_dsigdEtj1     ->fill(EtJet1);
        _h_DIS_dsigdMX       ->fill(sqrt(M2X)/GeV);
        _h_DIS_dsigdDeltaEta ->fill(deltaEtaJets);
        _h_DIS_dsigdAvgEta   ->fill(avgEtaJets);
      }

    }

    // Finalize
    void finalize() {
      // Normalise to cross section
      // Remember to manually scale the cross section afterwards with
      // the number of rejected events.
      const double norm = crossSection()/picobarn/sumOfWeights();

      scale(_h_PHO_sig_sqrts,     norm);
      scale(_h_PHO_dsigdz,        norm);
      scale(_h_PHO_dsigdxPom,     norm);
      scale(_h_PHO_dsigdy,        norm);
      scale(_h_PHO_dsigdxGam,     norm);
      scale(_h_PHO_dsigdEtj1,     norm);
      scale(_h_PHO_dsigdMX,       norm);
      scale(_h_PHO_dsigdDeltaEta, norm);
      scale(_h_PHO_dsigdAvgEta,   norm);

      scale(_h_DIS_sig_sqrts,     norm);
      scale(_h_DIS_dsigdz,        norm);
      scale(_h_DIS_dsigdxPom,     norm);
      scale(_h_DIS_dsigdy,        norm);
      scale(_h_DIS_dsigdQ2,       norm);
      scale(_h_DIS_dsigdEtj1,     norm);
      scale(_h_DIS_dsigdMX,       norm);
      scale(_h_DIS_dsigdDeltaEta, norm);
      scale(_h_DIS_dsigdAvgEta,   norm);

      if (_h_DIS_sig_sqrts->numEntries() != 0)
        divide(_h_PHO_sig_sqrts, _h_DIS_sig_sqrts, _h_PHODIS_sqrts);
      if (_h_DIS_dsigdDeltaEta->numEntries() != 0)
        divide(_h_PHO_dsigdDeltaEta, _h_DIS_dsigdDeltaEta, _h_PHODIS_deltaEta);
      if (_h_DIS_dsigdy->numEntries() != 0)
        divide(_h_PHO_dsigdy, _h_DIS_dsigdy, _h_PHODIS_y);
      if (_h_DIS_dsigdz->numEntries() != 0)
        divide(_h_PHO_dsigdz, _h_DIS_dsigdz, _h_PHODIS_z);
      if (_h_DIS_dsigdEtj1->numEntries() != 0)
        divide(_h_PHO_dsigdEtj1, _h_DIS_dsigdEtj1, _h_PHODIS_Etj1);

      const double dPHO = numEvents();
      MSG_INFO("H1_2015_I1343110");
      MSG_INFO("Cross section = " << crossSection()/picobarn << " pb");
      MSG_INFO("Number of events = " << numEvents() << ", sumW = " << sumOfWeights());
      MSG_INFO("Number of PHO = " << nPHO << ", number of DIS = " << nDIS);
      MSG_INFO("Events passing event setup     = " << nVeto0 << " (" << nVeto0/dPHO * 100. << "%)" );
      MSG_INFO("Events passing electron veto   = " << nVeto1 << " (" << nVeto1/dPHO * 100. << "%)" );
      MSG_INFO("Events passing t veto          = " << nVeto2 << " (" << nVeto2/dPHO * 100. << "%)" );
      MSG_INFO("Events passing xPom            = " << nVeto3 << " (" << nVeto3/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet Et   veto   = " << nVeto4 << " (" << nVeto4/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet eta veto    = " << nVeto5 << " (" << nVeto5/dPHO * 100. << "%)" );
      MSG_INFO("Events passing zPom veto       = " << nVeto6 << " (" << nVeto6/dPHO * 100. << "%)" );

    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    // Book histograms from REF data
    BinnedHistoPtr<int> _h_PHO_sig_sqrts;
    BinnedHistoPtr<int> _h_DIS_sig_sqrts;
    BinnedEstimatePtr<int> _h_PHODIS_sqrts;

    Histo1DPtr _h_DIS_dsigdz;
    Histo1DPtr _h_DIS_dsigdxPom;
    Histo1DPtr _h_DIS_dsigdy;
    Histo1DPtr _h_DIS_dsigdQ2;
    Histo1DPtr _h_DIS_dsigdEtj1;
    Histo1DPtr _h_DIS_dsigdMX;
    Histo1DPtr _h_DIS_dsigdDeltaEta;
    Histo1DPtr _h_DIS_dsigdAvgEta;

    Histo1DPtr _h_PHO_dsigdz;
    Histo1DPtr _h_PHO_dsigdxPom;
    Histo1DPtr _h_PHO_dsigdy;
    Histo1DPtr _h_PHO_dsigdxGam;
    Histo1DPtr _h_PHO_dsigdEtj1;
    Histo1DPtr _h_PHO_dsigdMX;
    Histo1DPtr _h_PHO_dsigdDeltaEta;
    Histo1DPtr _h_PHO_dsigdAvgEta;

    Estimate1DPtr _h_PHODIS_deltaEta;
    Estimate1DPtr _h_PHODIS_y;
    Estimate1DPtr _h_PHODIS_z;
    Estimate1DPtr _h_PHODIS_Etj1;
    /// @}

    bool isPHO;
    int  nVeto0, nVeto1, nVeto2, nVeto3, nVeto4, nVeto5, nVeto6;
    int  nPHO, nDIS;
  };

  RIVET_DECLARE_PLUGIN(H1_2015_I1343110);

}
