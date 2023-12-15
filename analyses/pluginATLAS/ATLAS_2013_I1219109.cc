// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {


  /// @brief ATLAS W+b measurement
  class ATLAS_2013_I1219109: public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1219109);


    void init() {

      // Get options from the option system
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;


      // Initialise and register projections
      declare("MET", MissingMomentum());
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
      LeptonFinder ef(0.1, cuts && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");
      LeptonFinder mf(0.1, cuts && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      // Jets
      VetoedFinalState jet_fs;
      jet_fs.addVetoOnThisFinalState(ef);
      jet_fs.addVetoOnThisFinalState(mf);
      FastJets fj(jet_fs, JetAlg::ANTIKT, 0.4);
      fj.useInvisibles();
      declare(fj, "Jets");
      declare(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "BHadrons");

      // Book histograms
      book(_njet     ,1, 1, 1); // dSigma / dNjet
      book(_jet1_bPt ,3, 1, 1); // dSigma / dBjetPt for Njet = 1
      book(_jet2_bPt ,8, 1, 1); // dSigma / dBjetPt for Njet = 2
    }


    void analyze(const Event& event) {

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Retrieve W boson candidate
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int iefound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV);
      const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
      const int imfound = closestMatchIndex(mus, pmiss, Kin::mass, 80.4*GeV);

      // Event selection and identify W lepton
      size_t nWel = (iefound >= 0);
      size_t nWmu = (imfound >= 0);
      if (_mode == 0 && !((nWmu == 1 && !nWel) || (!nWmu && nWel == 1)))  vetoEvent; // one W->munu OR W->elnu candidate, otherwise veto
      if (_mode == 1 && !(!nWmu && nWel == 1))  vetoEvent; // one W->elnu candidate, otherwise veto
      if (_mode == 2 && !(nWmu == 1 && !nWel))  vetoEvent; // one W->munu candidate, otherwise veto
      const Particle& lepton = (nWel > 0) ? es[iefound] : mus[imfound];

      // mT cut
      /// @todo Shouldn't this go earlier, in the lepton filtering?
      if (mT(lepton, pmiss) < 60*GeV) vetoEvent;

      // Count good jets, check if good jet contains B hadron
      const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();
      const Jets& jets = apply<JetFinder>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV);
      int goodjets = 0, bjets = 0;
      double bPt = 0.;
      for (const Jet& j : jets) {
        if ( (j.abseta() < 2.1) && (deltaR(lepton, j) > 0.5) ) {
          // this jet passes the selection!
          ++goodjets;
          // Perform manual dR matching
          for (const Particle& b : bHadrons) {
            if ( deltaR(j, b) < 0.3 ) {
              // jet matched to B hadron!
              if (!bPt) bPt = j.pT()/GeV; // leading b-jet pT
              ++bjets; // count number of b-jets
              break;
            }
          }
        }
      }
      if (goodjets > 2) vetoEvent; // at most two jets
      if (!bjets) vetoEvent; // at least one of them b-tagged

      const double ncomb = 3.0; //< ?
      _njet->fill(goodjets);
      _njet->fill(ncomb);

      if (goodjets == 1) _jet1_bPt->fill(bPt);
      else if (goodjets == 2) _jet2_bPt->fill(bPt);
    }


    void finalize() {
      const double sf = _mode? 1.0 : 0.5;
      const double xs_pb = sf * crossSection() / picobarn  / sumOfWeights();
      const double xs_fb = sf * crossSection() / femtobarn / sumOfWeights();
      scale(_njet,     xs_pb);
      scale(_jet1_bPt, xs_fb);
      scale(_jet2_bPt, xs_fb);
    }


  private:

    size_t _mode;

    Histo1DPtr _njet;
    Histo1DPtr _jet1_bPt;
    Histo1DPtr _jet2_bPt;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1219109);

}
