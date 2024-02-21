// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// An MC analysis for studying the content of {W,Z} H(->bb) events
  class MC_VH2BB : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_VH2BB);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      _jetptcut = getOption<double>("PTJMIN", 30.0)*GeV;

      // Boson projections
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      DileptonFinder zeefinder(91.2*GeV, 0.2, cut && Cuts::abspid == PID::ELECTRON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zeefinder, "ZeeFinder");
      DileptonFinder zmmfinder(91.2*GeV, 0.2, cut && Cuts::abspid == PID::MUON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zmmfinder, "ZmmFinder");
      //
      LeptonFinder ef(cut && Cuts::abspid == PID::ELECTRON, 0.2);
      declare(ef, "Elecs");
      LeptonFinder mf(cut && Cuts::abspid == PID::MUON, 0.2);
      declare(ef, "Muons");
      declare(MissingMomentum(), "MET");

      // Jet projections
      VetoedFinalState jfs;
      jfs
      .addVetoOnThisFinalState(zeefinder)
      .addVetoOnThisFinalState(zmmfinder)
      .addVetoOnThisFinalState(ef)
      .addVetoOnThisFinalState(mf);
      declare(FastJets(jfs, JetAlg::ANTIKT, 0.4), "AntiKT04");
      declare(FastJets(jfs, JetAlg::ANTIKT, 0.5), "AntiKT05");
      declare(FastJets(jfs, JetAlg::ANTIKT, 0.6), "AntiKT06");

      /// Book histograms
      book(_h_jet_bb_Delta_eta ,"jet_bb_Delta_eta", 50, 0, 4);
      book(_h_jet_bb_Delta_phi ,"jet_bb_Delta_phi", 50, 0, 1);
      book(_h_jet_bb_Delta_pT ,"jet_bb_Delta_pT", 50,0, 500);
      book(_h_jet_bb_Delta_R ,"jet_bb_Delta_R", 50, 0, 5);
      book(_h_jet_b_jet_eta ,"jet_b_jet_eta", 50, -4, 4);
      book(_h_jet_b_jet_multiplicity ,"jet_b_jet_multiplicity", 11, -0.5, 10.5);
      book(_h_jet_b_jet_phi ,"jet_b_jet_phi", 50, 0, 1);
      book(_h_jet_b_jet_pT ,"jet_b_jet_pT", 50, 0, 500);
      book(_h_jet_H_eta_using_bb ,"jet_H_eta_using_bb", 50, -4, 4);
      book(_h_jet_H_mass_using_bb ,"jet_H_mass_using_bb", 50, 50, 200);
      book(_h_jet_H_phi_using_bb ,"jet_H_phi_using_bb", 50, 0, 1);
      book(_h_jet_H_pT_using_bb ,"jet_H_pT_using_bb", 50, 0, 500);
      book(_h_jet_eta ,"jet_eta", 50, -4, 4);
      book(_h_jet_multiplicity ,"jet_multiplicity", 11, -0.5, 10.5);
      book(_h_jet_phi ,"jet_phi", 50, 0, 1);
      book(_h_jet_pT ,"jet_pT", 50, 0, 500);
      book(_h_jet_VBbb_Delta_eta ,"jet_VBbb_Delta_eta", 50, 0, 4);
      book(_h_jet_VBbb_Delta_phi ,"jet_VBbb_Delta_phi", 50, 0, 1);
      book(_h_jet_VBbb_Delta_pT ,"jet_VBbb_Delta_pT", 50, 0, 500);
      book(_h_jet_VBbb_Delta_R ,"jet_VBbb_Delta_R", 50, 0, 8);

      book(_h_VB_eta ,"VB_eta", 50, -4, 4);
      book(_h_VB_mass ,"VB_mass", 50, 60, 110);
      book(_h_Z_multiplicity ,"Z_multiplicity", 11, -0.5, 10.5);
      book(_h_W_multiplicity ,"W_multiplicity", 11, -0.5, 10.5);
      book(_h_VB_phi ,"VB_phi", 50, 0, 1);
      book(_h_VB_pT ,"VB_pT", 50, 0, 500);

      book(_h_jet_bVB_angle_Hframe ,"jet_bVB_angle_Hframe", 50, 0, 1);
      book(_h_jet_bb_angle_Hframe ,"jet_bb_angle_Hframe", 50, 0, 1);
      book(_h_jet_bVB_cosangle_Hframe ,"jet_bVB_cosangle_Hframe", 50, -1, 1);
      book(_h_jet_bb_cosangle_Hframe ,"jet_bb_cosangle_Hframe", 50, -1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get vector bosons -- only one is expected
      const DileptonFinder& zeefinder = apply<DileptonFinder>(event, "ZeeFinder");
      const DileptonFinder& zmmfinder = apply<DileptonFinder>(event, "ZmmFinder");
      Particles vectorBosons = zeefinder.bosons() + zmmfinder.bosons();
      const int numZ = vectorBosons.size();
      // Now the W's...
      const P4& pmiss = apply<MissingMomentum>(event, "MET").missingMom();
      if (pmiss.pT() > 25*GeV) {
        const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
        const int iefound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
        if (iefound >= 0) {
          const Particle we(PID::WBOSON * sign(es[iefound].charge()), es[iefound].mom() + pmiss);
          vectorBosons.push_back(we);
        }
        //
        const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
        const int imfound = closestMatchIndex(mus, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
        if (imfound >= 0) {
          const Particle wm(PID::WBOSON * sign(mus[imfound].charge()), mus[imfound].mom() + pmiss);
          vectorBosons.push_back(wm);
        }
      }
      const int numW = vectorBosons.size() - numZ;
      _h_Z_multiplicity->fill(numZ);
      _h_W_multiplicity->fill(numW);

      const Jets jets = apply<FastJets>(event, "AntiKT04").jetsByPt(Cuts::pT > _jetptcut);
      _h_jet_multiplicity->fill(jets.size());

      // Identify the b-jets
      Jets bjets;
      for (const Jet& jet : jets) {
        const double jetEta = jet.eta();
        const double jetPhi = jet.phi();
        const double jetPt = jet.pT();
        _h_jet_eta->fill(jetEta);
        _h_jet_phi->fill(jetPhi/2/M_PI);
        _h_jet_pT->fill(jetPt/GeV);

        if (jet.bTagged() && jet.pT() > _jetptcut) {
          bjets.push_back(jet);
          _h_jet_b_jet_eta->fill(jetEta);
          _h_jet_b_jet_phi->fill(jetPhi/2/M_PI);
          _h_jet_b_jet_pT->fill(jetPt);
        }
      }
      _h_jet_b_jet_multiplicity->fill(bjets.size());

      // Plot vector-boson properties
      for (const Particle& v : vectorBosons) {
        _h_VB_phi->fill(v.phi()/2/M_PI);
        _h_VB_pT->fill(v.pT());
        _h_VB_eta->fill(v.eta());
        _h_VB_mass->fill(v.mass());
      }

      // The rest of the analysis requires at least 1 b-jet
      if (bjets.empty()) vetoEvent;

      // Construct Higgs candidates from pairs of b-jets
      for (size_t i = 0; i < bjets.size()-1; ++i) {
        for (size_t j = i+1; j < bjets.size(); ++j) {
          const Jet& jet1 = bjets[i];
          const Jet& jet2 = bjets[j];

          const double deltaEtaJJ = fabs(jet1.eta() - jet2.eta());
          const double deltaPhiJJ = deltaPhi(jet1.momentum(), jet2.momentum());
          const double deltaRJJ = deltaR(jet1.momentum(), jet2.momentum());
          const double deltaPtJJ = fabs(jet1.pT() - jet2.pT());
          _h_jet_bb_Delta_eta->fill(deltaEtaJJ);
          _h_jet_bb_Delta_phi->fill(deltaPhiJJ/M_PI);
          _h_jet_bb_Delta_pT->fill(deltaPtJJ);
          _h_jet_bb_Delta_R->fill(deltaRJJ);

          const FourMomentum phiggs = jet1.momentum() + jet2.momentum();
          _h_jet_H_eta_using_bb->fill(phiggs.eta());
          _h_jet_H_mass_using_bb->fill(phiggs.mass());
          _h_jet_H_phi_using_bb->fill(phiggs.phi()/2/M_PI);
          _h_jet_H_pT_using_bb->fill(phiggs.pT());

          for (const Particle& v : vectorBosons) {
            const double deltaEtaVH = fabs(phiggs.eta() - v.eta());
            const double deltaPhiVH = deltaPhi(phiggs, v.momentum());
            const double deltaRVH = deltaR(phiggs, v.momentum());
            const double deltaPtVH = fabs(phiggs.pT() - v.pT());
            _h_jet_VBbb_Delta_eta->fill(deltaEtaVH);
            _h_jet_VBbb_Delta_phi->fill(deltaPhiVH/M_PI);
            _h_jet_VBbb_Delta_pT->fill(deltaPtVH);
            _h_jet_VBbb_Delta_R->fill(deltaRVH);

            // Calculate boost angles
            const vector<double> angles = boostAngles(jet1.momentum(), jet2.momentum(), v.momentum());
            _h_jet_bVB_angle_Hframe->fill(angles[0]/M_PI);
            _h_jet_bb_angle_Hframe->fill(angles[1]/M_PI);
            _h_jet_bVB_cosangle_Hframe->fill(cos(angles[0]));
            _h_jet_bb_cosangle_Hframe->fill(cos(angles[1]));
          }

        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_jet_bb_Delta_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_bb_Delta_phi, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_bb_Delta_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_bb_Delta_R, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_b_jet_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_b_jet_multiplicity, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_b_jet_phi, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_b_jet_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_H_eta_using_bb, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_H_mass_using_bb, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_H_phi_using_bb, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_H_pT_using_bb, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_multiplicity, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_phi, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_VBbb_Delta_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_VBbb_Delta_phi, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_VBbb_Delta_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_VBbb_Delta_R, crossSection()/picobarn/sumOfWeights());

      scale(_h_VB_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_VB_mass, crossSection()/picobarn/sumOfWeights());
      scale(_h_Z_multiplicity, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_multiplicity, crossSection()/picobarn/sumOfWeights());
      scale(_h_VB_phi, crossSection()/picobarn/sumOfWeights());
      scale(_h_VB_pT, crossSection()/picobarn/sumOfWeights());

      scale(_h_jet_bVB_angle_Hframe, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_bb_angle_Hframe, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_bVB_cosangle_Hframe, crossSection()/picobarn/sumOfWeights());
      scale(_h_jet_bb_cosangle_Hframe, crossSection()/picobarn/sumOfWeights());
    }


    /// This should take in the four-momenta of two b's (jets/hadrons) and a vector boson, for the process VB*->VBH with H->bb
    /// It should return the smallest angle between the virtual vector boson and one of the b's, in the rest frame of the Higgs boson.
    /// It should also return (as the second element of the vector) the angle between the b's, in the rest frame of the Higgs boson.
    vector<double> boostAngles(const FourMomentum& b1, const FourMomentum& b2, const FourMomentum& vb) {
      const FourMomentum higgsMomentum = b1 + b2;
      const FourMomentum virtualVBMomentum = higgsMomentum + vb;
      const LorentzTransform lt = LorentzTransform::mkFrameTransformFromBeta(higgsMomentum.betaVec());

      const FourMomentum virtualVBMomentumBOOSTED = lt.transform(virtualVBMomentum);
      const FourMomentum b1BOOSTED = lt.transform(b1);
      const FourMomentum b2BOOSTED = lt.transform(b2);

      const double angle1 = b1BOOSTED.angle(virtualVBMomentumBOOSTED);
      const double angle2 = b2BOOSTED.angle(virtualVBMomentumBOOSTED);

      const double anglebb = mapAngle0ToPi(b1BOOSTED.angle(b2BOOSTED));

      vector<double> rtn;
      rtn.push_back(angle1 < angle2 ? angle1 : angle2);
      rtn.push_back(anglebb);
      return rtn;
    }

    /// @}


  private:

    // Jet pT cut
    double _jetptcut = 20*GeV;

    /// @name Histograms
    /// @{
    Histo1DPtr _h_Z_multiplicity, _h_W_multiplicity;
    Histo1DPtr _h_jet_bb_Delta_eta, _h_jet_bb_Delta_phi, _h_jet_bb_Delta_pT, _h_jet_bb_Delta_R;
    Histo1DPtr _h_jet_b_jet_eta, _h_jet_b_jet_multiplicity, _h_jet_b_jet_phi, _h_jet_b_jet_pT;
    Histo1DPtr _h_jet_H_eta_using_bb, _h_jet_H_mass_using_bb, _h_jet_H_phi_using_bb, _h_jet_H_pT_using_bb;
    Histo1DPtr _h_jet_eta, _h_jet_multiplicity, _h_jet_phi, _h_jet_pT;
    Histo1DPtr _h_jet_VBbb_Delta_eta, _h_jet_VBbb_Delta_phi, _h_jet_VBbb_Delta_pT, _h_jet_VBbb_Delta_R;
    Histo1DPtr _h_VB_eta, _h_VB_mass, _h_VB_phi, _h_VB_pT;
    Histo1DPtr _h_jet_bVB_angle_Hframe, _h_jet_bb_angle_Hframe, _h_jet_bVB_cosangle_Hframe, _h_jet_bb_cosangle_Hframe;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_VH2BB);

}
