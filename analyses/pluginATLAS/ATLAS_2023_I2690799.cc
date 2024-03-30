// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Measurement of electroweak 4 lepton + 2 jet production
  class ATLAS_2023_I2690799 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2690799);

    /// A Z Boson structrue:
    struct Z_can {
       int lep_1, lep_2;
       double dZ;
       FourMomentum p;
    };

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections


      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON);
      PromptFinalState muons(Cuts::abspid == PID::MUON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.05,
      // and apply some fiducial cuts on the dressed leptons
      Cut cuts_el = (Cuts::pT > 7*GeV) && ( Cuts::abseta < 2.47);
      Cut cuts_mu = (Cuts::pT > 5*GeV) && ( Cuts::abseta < 2.7);

      LeptonFinder dressed_electrons(electrons, photons, 0.1, cuts_el);
      declare(dressed_electrons, "DressedElectrons");

      LeptonFinder dressed_muons(muons, photons, 0.1, cuts_mu);
      declare(dressed_muons, "DressedMuons");

      // We define AntiKt4TruthWZJets here, the leptons here are not used for next step selection.
      // Muons
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT); // true = use muons from prompt tau decays
      LeptonFinder all_dressed_mu(bare_mu, photons, 0.1, Cuts::abseta < 2.5);
      //
      // Electrons
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT); // true = use electrons from prompt tau decays
      LeptonFinder all_dressed_el(bare_el, photons, 0.1, Cuts::abseta < 2.5);

      VetoedFinalState vfs(FinalState(Cuts::abseta < 4.5));
      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // AntiKt4TruthWZJets
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);
      FastJets jetfs(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jetfs, "jets");

      // Book histograms
      book( _h["m_jj_SR"],     7, 1, 1);
      book( _h["m_jj_CR"],     8, 1, 1);
      book( _h["m_4l_SR"],     9, 1, 1);
      book( _h["m_4l_CR"],    10, 1, 1);
      book( _h["pt_4l_SR"],   11, 1, 1);
      book( _h["pt_4l_CR"],   12, 1, 1);
      book( _h["dphi_SR"],    13, 1, 1);
      book( _h["dphi_CR"],    14, 1, 1);
      book( _h["dY_SR"],      15, 1, 1);
      book( _h["dY_CR"],      16, 1, 1);
      book( _h["costs1_SR"],  17, 1, 1);
      book( _h["costs1_CR"],  18, 1, 1);
      book( _h["costs3_SR"],  19, 1, 1);
      book( _h["costs3_CR"],  20, 1, 1);
      book( _h["pt_jj_SR"],   21, 1, 1);
      book( _h["pt_jj_CR"],   22, 1, 1);
      book( _h["pt_4ljj_SR"], 23, 1, 1);
      book( _h["pt_4ljj_CR"], 24, 1, 1);
      book( _h["st_4ljj_SR"], 25, 1, 1);
      book( _h["st_4ljj_CR"], 26, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      DressedLeptons muons = apply<LeptonFinder>(event, "DressedMuons").dressedLeptons();
      DressedLeptons elecs = apply<LeptonFinder>(event, "DressedElectrons").dressedLeptons();
      DressedLeptons leps;
      leps.insert( leps.end(), muons.begin(), muons.end() );
      leps.insert( leps.end(), elecs.begin(), elecs.end() );

      size_t n_leps = leps.size();
      if (n_leps < 4) vetoEvent;

      //leading and sub-leading leps pt
      sort(leps.begin(), leps.end(), cmpMomByPt);
      if (leps[0].pT() < 20*GeV || leps[1].pT() < 20*GeV) vetoEvent;
      //lepton seperation
      double dR_min = 9999;
      for (size_t i = 0; i < n_leps - 1; ++i) {
        for (size_t j = i + 1; j < n_leps; ++j) {
          double dR = deltaR(leps[i], leps[j]);
          if (dR < dR_min)  dR_min = dR;
        }
      }
      if (dR_min < 0.05) vetoEvent;


      //ZZ pairing
      size_t sfos_pair = 0;
      double minm2l    = 9999.;
      vector<Z_can> z;
      Particle *l1, *l3;
      for (size_t i = 0; i < n_leps - 1; ++i) {
        for (size_t j = i + 1; j < n_leps; ++j) {

          if (leps[i].pid() != -1*leps[j].pid())  continue;
          ++sfos_pair;
          FourMomentum twoL = leps[i].momentum() + leps[j].momentum();
          if (twoL.mass() < minm2l) minm2l = twoL.mass();

          double deltaM = fabs( twoL.mass() - Zmass_PDG);
          Z_can tmp_z;
          tmp_z.dZ = deltaM;
          tmp_z.p  = twoL;
          if (leps[i].pid() < 0) {
            tmp_z.lep_1 = i;
            tmp_z.lep_2 = j;
          } else{
            tmp_z.lep_1 = j;
            tmp_z.lep_2 = i;
          }
          z.push_back(tmp_z);
        }
      }
      if (minm2l < 5*GeV) vetoEvent; //All SFOS ll pairs should have m2l > 5 GeV.
      if (sfos_pair < 2)  vetoEvent;


      std::sort(z.begin(), z.end(), [](const Z_can& z1, const Z_can& z2) { return z1.dZ < z2.dZ;});

      int z1_id = -1;
      int z2_id = -1;
      FourMomentum fourL;
      bool found_4l = false;
      for (size_t i = 0; i < z.size() - 1; ++i) {
        for (size_t j = i + 1; j < z.size(); ++j) {
          if (z[i].lep_1 == z[j].lep_1 ||
              z[i].lep_1 == z[j].lep_2 ||
              z[i].lep_2 == z[j].lep_1 ||
              z[i].lep_2 == z[j].lep_2) continue;
          fourL = z[i].p + z[j].p;
          if (fourL.mass() > 130*GeV) {
            z1_id = i;
            z2_id = j;
            found_4l =true;
            break;
          }
        }
        if (found_4l)  break;
      }
      if (!found_4l) vetoEvent;

      // To judge which is the leading pair

      FourMomentum twoL_1, twoL_2;

      if (abs(z[z1_id].p.rap()) > abs(z[z2_id].p.rap())){
        twoL_1 = z[z1_id].p;
        twoL_2 = z[z2_id].p;
        int i = z[z1_id].lep_1;
        int j = z[z2_id].lep_1;
        l1 = &leps[i];
        l3 = &leps[j];
      } else {
        twoL_1 = z[z2_id].p;
        twoL_2 = z[z1_id].p;
        int i = z[z2_id].lep_1;
        int j = z[z1_id].lep_1;
        l1 = &leps[i];
        l3 = &leps[j];
      }

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.5);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].pt() < 40*GeV) vetoEvent;
      FourMomentum twoJ;
      double dY = -1;
      int jet_id1 = -1, jet_id2 = -1;
      bool found_2j = false;
      for (size_t i = 0; i < jets.size() - 1; ++i) {
        for (size_t j = i + 1; j < jets.size(); ++j) {
          if (jets[i].rap()*jets[j].rap() > 0)  continue; // First pair in A/C
          dY = fabs(jets[i].rap() - jets[j].rap());
          twoJ = jets[i].momentum() + jets[j].momentum();
          jet_id1 = i;
          jet_id2 = j;
          found_2j = true;
          break;
        }
        if (found_2j)  break;
      }
      if (!found_2j)  vetoEvent;
      if (dY < 2)     vetoEvent;
      if (twoJ.mass() < 300*GeV) vetoEvent;

      // Calculate Observables!

      double phi1, phi2;
      if (jets[jet_id1].rap() > jets[jet_id2].rap()) {
        phi1 = jets[jet_id1].phi(); phi2 = jets[jet_id2].phi();
      }
      else {
        phi1 = jets[jet_id2].phi(); phi2 = jets[jet_id1].phi();
      }
      const double dphi = mapAngleMPiToPi(phi1 - phi2);
      const double mjj  = twoJ.mass();
      const double m4l  = fourL.mass();
      const double pt4l = fourL.pT();
      const double ptjj = twoJ.pT();
      const double pt4ljj = (fourL + twoJ).pT();
      const double st4ljj = twoL_1.pT() + twoL_2.pT() + jets[jet_id1].pT() + jets[jet_id2].pT();
      const double centrality =   abs(fourL.rap() - 0.5*(jets[jet_id1].rap() + jets[jet_id2].rap())) \
                                / abs(jets[jet_id1].rap() - jets[jet_id2].rap());

      //negetive leps1 in Z rest frame work
      Vector3 beta2lcom_1 = twoL_1.betaVec();
      LorentzTransform com2lboost_1 = LorentzTransform::mkFrameTransformFromBeta(beta2lcom_1);
      FourMomentum l1_2lcom = com2lboost_1.transform(l1->momentum());
      const double costs1 = cos(twoL_1.angle(l1_2lcom));


      //negetive leps3 in Z rest frame work
      Vector3 beta2lcom_2 = twoL_2.betaVec();
      LorentzTransform com2lboost_2 = LorentzTransform::mkFrameTransformFromBeta(beta2lcom_2);
      FourMomentum l3_2lcom = com2lboost_2.transform(l3->momentum());
      const double costs3 = cos(twoL_2.angle(l3_2lcom));

      const string suff = (centrality < 0.4)? "_SR" : "_CR";
      _h["dphi"+suff]->fill(dphi);
      _h["m_jj"+suff]->fill(mjj/GeV);
      _h["dY"+suff]->fill(dY);
      _h["m_4l"+suff]->fill(m4l/GeV);
      _h["pt_4l"+suff]->fill(pt4l/GeV);
      _h["pt_jj"+suff]->fill(ptjj/GeV);
      _h["pt_4ljj"+suff]->fill(pt4ljj/GeV);
      _h["st_4ljj"+suff]->fill(st4ljj/GeV);
      _h["costs1"+suff]->fill(costs1/GeV);
      _h["costs3"+suff]->fill(costs3/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = crossSection()/sumOfWeights()/femtobarn;
      scale(_h, norm);

    }

    private:

      map<string, Histo1DPtr> _h;
      const double Zmass_PDG = 91.1876;
  };


  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2690799);

}
