// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class ATLAS_2011_I944826 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I944826);

    /// Book histograms and initialise projections before the run
    void init() {

      UnstableParticles ufs(Cuts::pT > 100*MeV);
      declare(ufs, "UFS");

      ChargedFinalState  mbts(Cuts::absetaIn(2.09, 3.84));
      declare(mbts, "MBTS");

      IdentifiedFinalState nstable(Cuts::abseta < 2.5 && Cuts::pT >= 100*MeV);
      nstable.acceptIdPair(PID::ELECTRON)
             .acceptIdPair(PID::MUON)
             .acceptIdPair(PID::PIPLUS)
             .acceptIdPair(PID::KPLUS)
             .acceptIdPair(PID::PROTON);
      declare(nstable, "nstable");


      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        bool is900(en == "900");
        size_t offset = is900? 3 : 0;
        book(_h[en+"Ks_pT"],       1+offset, 1, 1);
        book(_h[en+"Ks_y"],        2+offset, 1, 1);
        book(_h[en+"Ks_mult"],     3+offset, 1, 1);
        book(_h[en+"L_pT"],        7+offset, 1, 1);
        book(_h[en+"L_y"],         8+offset, 1, 1);
        book(_h[en+"L_mult"],      9+offset, 1, 1);
        if (is900)  offset = 2;
        book(_e[en+"v_y"],  13+offset, 1, 1);
        book(_e[en+"v_pT"], 14+offset, 1, 1);
        //
        book(_h[en+"lambda_v_y"],     "TMP/lambda_v_y_"+en,     is900? 5 : 10, 0.0, 2.5);
        book(_h[en+"lambdabar_v_y"],  "TMP/lambdabar_v_y_"+en,  is900? 5 : 10, 0.0, 2.5);
        book(_h[en+"lambda_v_pT"],    "TMP/lambda_v_pT_"+en,    is900? 8 : 18, 0.5, is900? 3.7 : 4.1);
        book(_h[en+"lambdabar_v_pT"], "TMP/lambdabar_v_pT_"+en, is900? 8 : 18, 0.5, is900? 3.7 : 4.1);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    // This function is required to impose the flight time cuts on Kaons and Lambdas
    double getPerpFlightDistance(const Rivet::Particle& p) {
      ConstGenParticlePtr genp = p.genParticle();
      ConstGenVertexPtr prodV = genp->production_vertex();
      ConstGenVertexPtr decV  = genp->end_vertex();
      RivetHepMC::FourVector prodPos = prodV->position();
      if (decV) {
        const RivetHepMC::FourVector decPos = decV->position();
        double dy = prodPos.y() - decPos.y();
        double dx = prodPos.x() - decPos.x();
        return add_quad(dx, dy);
      }
      return numeric_limits<double>::max();
    }


    bool daughtersSurviveCuts(const Rivet::Particle& p) {
      // We require the Kshort or Lambda to decay into two charged
      // particles with at least pT = 100 MeV inside acceptance region
      ConstGenParticlePtr genp = p.genParticle();
      ConstGenVertexPtr decV  = genp->end_vertex();
      bool decision = true;

      if (!decV) return false;
      if (HepMCUtils::particles(decV, Relatives::CHILDREN).size() == 2) {
        std::vector<double> pTs;
        std::vector<int> charges;
        std::vector<double> etas;
        for(ConstGenParticlePtr gp: HepMCUtils::particles(decV, Relatives::CHILDREN)) {
          pTs.push_back(gp->momentum().perp());
          etas.push_back(fabs(gp->momentum().eta()));
          charges.push_back( Rivet::PID::charge3(gp->pdg_id()) );
          // gp->print();
        }
        if ( (pTs[0]/Rivet::GeV < 0.1) || (pTs[1]/Rivet::GeV < 0.1) ) {
          decision = false;
          MSG_DEBUG("Failed pT cut: " << pTs[0]/Rivet::GeV << " " << pTs[1]/Rivet::GeV);
        }
        if ( etas[0] > 2.5 || etas[1] > 2.5 ) {
          decision = false;
          MSG_DEBUG("Failed eta cut: " << etas[0] << " " << etas[1]);
        }
        if ( charges[0] * charges[1] >= 0 ) {
          decision = false;
          MSG_DEBUG("Failed opposite charge cut: " << charges[0] << " " << charges[1]);
        }
      }
      else {
        decision = false;
        MSG_DEBUG("Failed nDaughters cut: " << HepMCUtils::particles(decV, Relatives::CHILDREN).size());
      }

      return decision;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // ATLAS MBTS trigger requirement of at least one hit in either hemisphere
      if (apply<FinalState>(event, "MBTS").size() < 1) {
        MSG_DEBUG("Failed trigger cut");
        vetoEvent;
      }

      // Veto event also when we find less than 2 particles in the acceptance region of type 211,2212,11,13,321
      if (apply<FinalState>(event, "nstable").size() < 2) {
        MSG_DEBUG("Failed stable particle cut");
        vetoEvent;
      }

      // This ufs holds all the Kaons and Lambdas
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      // Some conters
      int n_KS0 = 0, n_LAMBDA = 0;

      // Particle loop
      for (const Particle& p : ufs.particles()) {

        // General particle quantities
        const double pT = p.pT();
        const double y = p.rapidity();
        const PdgId apid = p.abspid();

        double flightd = 0.0;

        // Look for Kaons, Lambdas
        switch (apid) {

        case PID::K0S:
          flightd = getPerpFlightDistance(p);
          if (!inRange(flightd/mm, 4., 450.) ) {
            MSG_DEBUG("Kaon failed flight distance cut:" << flightd);
            break;
          }
          if (daughtersSurviveCuts(p) ) {
            _h[_sqs+"Ks_y"] ->fill(y);
            _h[_sqs+"Ks_pT"]->fill(pT/GeV);
            ++n_KS0;
          }
          break;

        case PID::LAMBDA:
          if (pT < 0.5*GeV) { // Lambdas have an additional pT cut of 500 MeV
            MSG_DEBUG("Lambda failed pT cut:" << pT/GeV << " GeV");
            break;
          }
          flightd = getPerpFlightDistance(p);
          if (!inRange(flightd/mm, 17., 450.)) {
            MSG_DEBUG("Lambda failed flight distance cut:" << flightd/mm << " mm");
            break;
          }
          if ( daughtersSurviveCuts(p) ) {
            if (p.pid() == PID::LAMBDA) {
              _h[_sqs+"lambda_v_y"]->fill(fabs(y));
              _h[_sqs+"lambda_v_pT"]->fill(pT/GeV);
              _h[_sqs+"L_y"]->fill(y);
              _h[_sqs+"L_pT"]->fill(pT/GeV);
              ++n_LAMBDA;
            }
            else if (p.pid() == -PID::LAMBDA) {
              _h[_sqs+"lambdabar_v_y"]->fill(fabs(y));
              _h[_sqs+"lambdabar_v_pT"]->fill(pT/GeV);
            }
          }
          break;
        }
      }

      // Fill multiplicity histos
      _h[_sqs+"Ks_mult"]->fill(n_KS0);
      _h[_sqs+"L_mult"]->fill(n_LAMBDA);
    }



    /// Normalise histograms etc., after the run
    void finalize() {
      // Division of histograms to obtain lambda_bar/lambda ratios
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));

        if (_h[en+"L_pT"]->sumW())  scale(_h[en+"L_y"], 1.0 / _h[en+"L_pT"]->sumW());
        divide(_h[en+"lambdabar_v_y"],  _h[en+"lambda_v_y"],  _e[en+"v_y"]);
        divide(_h[en+"lambdabar_v_pT"], _h[en+"lambda_v_pT"], _e[en+"v_pT"]);
      }
      normalize(_h);

    }


  private:

    /// @name Persistent histograms
    /// @{
    map<string,Histo1DPtr> _h;
    map<string,Estimate1DPtr> _e;

    string _sqs = "";
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2011_I944826);

}
