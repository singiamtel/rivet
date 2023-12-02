// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// Differential cross-sections for W boson and jets in proton-proton collisions at 13 TeV
  class CMS_2017_I1610623 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1610623);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      LeptonFinder lf(0.1, Cuts::abseta < 2.4 && Cuts::abspid == PID::MUON);
      declare(lf, "Leptons");
      declare(MissingMom(), "MET");

      // Define veto FS
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(lf);
      vfs.addVetoPairId(PID::MUON);
      vfs.vetoNeutrinos();

      FastJets fastjets(vfs, JetAlg::ANTIKT, 0.4);
      declare(fastjets, "Jets");


      book(_hist_Mult_exc      ,"d01-x01-y01");
      book(_hist_inc_WJetMult  ,"d02-x01-y01");

      book(_hist_JetPt1j ,"d03-x01-y01");
      book(_hist_JetPt2j ,"d04-x01-y01");
      book(_hist_JetPt3j ,"d05-x01-y01");
      book(_hist_JetPt4j ,"d06-x01-y01");

      book(_hist_JetRap1j ,"d07-x01-y01");
      book(_hist_JetRap2j ,"d08-x01-y01");
      book(_hist_JetRap3j ,"d09-x01-y01");
      book(_hist_JetRap4j ,"d10-x01-y01");

      book(_hist_Ht_1j ,"d11-x01-y01");
      book(_hist_Ht_2j ,"d12-x01-y01");
      book(_hist_Ht_3j ,"d13-x01-y01");
      book(_hist_Ht_4j ,"d14-x01-y01");

      book(_hist_dphij1mu_1j , "d15-x01-y01");
      book(_hist_dphij2mu_2j , "d16-x01-y01");
      book(_hist_dphij3mu_3j , "d17-x01-y01");
      book(_hist_dphij4mu_4j , "d18-x01-y01");

      book(_hist_dRmuj_1j , "d19-x01-y01");

    }


    /// Used for filling inc Njets histo
    void _fill(Histo1DPtr& _histJetMult, vector<FourMomentum>& finaljet_list) {
      _histJetMult->fill(0);
      for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
        if (i==6) break;
        _histJetMult->fill(i+1);  // inclusive multiplicity
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Identify the closest-matching mu+MET to m == mW
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      const Particles& mus = apply<LeptonFinder>(event, "Leptons").particles();
      const Particles mus_mtfilt = select(mus, [&](const Particle& m){ return mT(m, pmiss) > 50*GeV; });
      const int imfound = closestMatchIndex(mus_mtfilt, pmiss, Kin::mass, 80.4*GeV);

      // Make cuts on the identified mT and lepton
      if (imfound < 0) vetoEvent;
      const FourMomentum& lepton0 = mus_mtfilt[imfound].momentum();
      if (lepton0.abseta() > 2.4 || lepton0.pT() < 25.0*GeV) vetoEvent;

      // Obtain the jets
      vector<FourMomentum> finaljet_list, jet100_list;
      double HT = 0.0;
      const Jets& js = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 2.4);
      for (const Jet& j : js) {
        const double jpt = j.momentum().pT();
        if (deltaR(lepton0, j.momentum()) > 0.4) {
          if (jpt > 30.0*GeV) {
            finaljet_list.push_back(j.momentum());
            HT += j.momentum().pT();
          }
          if (jpt > 100.0*GeV) {
            jet100_list.push_back(j.momentum());
          }
        }
      }

      // Multiplicity plots
      _hist_Mult_exc->fill(finaljet_list.size());
      _fill(_hist_inc_WJetMult, finaljet_list);

      // dRmuj plot.
      double mindR(99999);
      if (jet100_list.size() >= 1) {
        for (unsigned ji = 0; ji < jet100_list.size(); ji++){
          double dr_(9999);
          dr_ = fabs(deltaR(lepton0, jet100_list[ji]));
          if (dr_ < mindR){
            mindR = dr_;
          }
        }
        if (jet100_list[0].pT() > 300.0*GeV) {
          _hist_dRmuj_1j->fill(mindR);
        }
      }

      if (finaljet_list.size()>=1) {
        _hist_JetPt1j->fill(finaljet_list[0].pT());
        _hist_JetRap1j->fill(fabs(finaljet_list[0].rap()));
        _hist_Ht_1j->fill(HT);
        _hist_dphij1mu_1j->fill(deltaPhi(finaljet_list[0].phi(), lepton0.phi()));
      }

      if (finaljet_list.size()>=2) {
        _hist_JetPt2j->fill(finaljet_list[1].pT());
        _hist_JetRap2j->fill(fabs(finaljet_list[1].rap()));
        _hist_Ht_2j->fill(HT);
        _hist_dphij2mu_2j->fill(deltaPhi(finaljet_list[1].phi(), lepton0.phi()));
      }

      if (finaljet_list.size()>=3) {
        _hist_JetPt3j->fill(finaljet_list[2].pT());
        _hist_JetRap3j->fill(fabs(finaljet_list[2].rap()));
        _hist_Ht_3j->fill(HT);
        _hist_dphij3mu_3j->fill(deltaPhi(finaljet_list[2].phi(), lepton0.phi()));
      }

      if (finaljet_list.size()>=4) {
        _hist_JetPt4j->fill(finaljet_list[3].pT());
        _hist_JetRap4j->fill(fabs(finaljet_list[3].rap()));
        _hist_Ht_4j->fill(HT);
        _hist_dphij4mu_4j->fill(deltaPhi(finaljet_list[3].phi(), lepton0.phi()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double crossec = !std::isnan(crossSectionPerEvent()) ? crossSection() : 61526.7*picobarn;
      if (std::isnan(crossSectionPerEvent())){
        MSG_INFO("No valid cross-section given, using NNLO xsec calculated by FEWZ " << crossec/picobarn << " pb");
      }

      scale(_hist_Mult_exc, crossec/picobarn/sumOfWeights());
      scale(_hist_inc_WJetMult, crossec/picobarn/sumOfWeights());

      scale(_hist_JetPt1j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetPt2j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetPt3j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetPt4j, crossec/picobarn/sumOfWeights());

      scale(_hist_JetRap1j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetRap2j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetRap3j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetRap4j, crossec/picobarn/sumOfWeights());

      scale(_hist_Ht_1j, crossec/picobarn/sumOfWeights());
      scale(_hist_Ht_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_Ht_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_Ht_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dphij1mu_1j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij2mu_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij3mu_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij4mu_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dRmuj_1j, crossec/picobarn/sumOfWeights());
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _hist_Mult_exc;
    Histo1DPtr _hist_inc_WJetMult;

    Histo1DPtr _hist_JetPt1j;
    Histo1DPtr _hist_JetPt2j;
    Histo1DPtr _hist_JetPt3j;
    Histo1DPtr _hist_JetPt4j;

    Histo1DPtr _hist_JetRap1j;
    Histo1DPtr _hist_JetRap2j;
    Histo1DPtr _hist_JetRap3j;
    Histo1DPtr _hist_JetRap4j;

    Histo1DPtr _hist_Ht_1j;
    Histo1DPtr _hist_Ht_2j;
    Histo1DPtr _hist_Ht_3j;
    Histo1DPtr _hist_Ht_4j;

    Histo1DPtr _hist_dphij1mu_1j;
    Histo1DPtr _hist_dphij2mu_2j;
    Histo1DPtr _hist_dphij3mu_3j;
    Histo1DPtr _hist_dphij4mu_4j;

    Histo1DPtr _hist_dRmuj_1j;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CMS_2017_I1610623);

}
