// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Average pT w.r.t. thrust and sphericity axes
  class PLUTO_1983_I191161 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1983_I191161);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      declare(ChargedFinalState(), "CFS");
      // Thrust
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      const Sphericity sphericity(fs);
      declare(sphericity, "Sphericity");

      sqs = 1.0;
      if (isCompatibleWithSqrtS(7.7))       sqs = "7.7";
      else if (isCompatibleWithSqrtS(9.4))  sqs = "9.4";
      else if (isCompatibleWithSqrtS(12.))  sqs = "12.0";
      else if (isCompatibleWithSqrtS(13.))  sqs = "13.0";
      else if (isCompatibleWithSqrtS(17.))  sqs = "17.0";
      else if (isCompatibleWithSqrtS(22.))  sqs = "22.0";
      else if (isCompatibleWithSqrtS(27.6)) sqs = "27.6";
      else if (isCompatibleWithSqrtS(30.8)) sqs = "30.0 - 31.6";
      else MSG_ERROR("Beam energy " << sqrtS() << " not supported!");
      book(_p_thrust_pt      , 1,1,1);
      book(_p_thrust_pt2     , 1,1,2);
      book(_p_thrust_sum_pt  , 1,1,3);
      book( _p_thrust_sum_pt2, 1,1,4);
      book(_p_sphere_pt      , 2,1,1);
      book( _p_sphere_pt2    , 2,1,2);
      book( _p_sphere_sum_pt , 2,1,3);
      book( _p_sphere_sum_pt2, 2,1,4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // at least 4 charged particles
      if (apply<ChargedFinalState>(event, "CFS").particles().size()<4) vetoEvent;
      // Sphericities
      MSG_DEBUG("Calculating sphericity");
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const FinalState & fs = apply<FinalState>(event, "FS");
      // remove pi0->gammagamma decay products and replace with the pi0s
      // needed to get the average pTs right
      Particles fsParticles;
      set<ConstGenParticlePtr> pi0;
      for(const Particle & p : fs.particles()) {
        if ((p.pid()!=PID::PHOTON && p.abspid()!=PID::ELECTRON)|| p.parents().empty() || p.parents()[0].pid()!=PID::PI0)
          fsParticles.push_back(p);
        else {
          if (pi0.find(p.parents()[0].genParticle())==pi0.end())
            fsParticles.push_back(p.parents()[0]);
        }
      }
      double pT_T_sum(0.),pT2_T_sum(0.);
      double pT_S_sum(0.),pT2_S_sum(0.);
      for(const Particle & p : fsParticles) {
        const Vector3 mom3 = p.p3();
        const double pTinT = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
        const double pT2_T = sqr(pTinT) + sqr(pToutT);
        const double pT2_S = sqr(pTinS) + sqr(pToutS);
        const double pT_T  = sqrt(pT2_T);
        const double pT_S  = sqrt(pT2_S);
        pT_T_sum  += sqrt(pT2_T);
        pT2_T_sum +=      pT2_T ;
        pT_S_sum  += sqrt(pT2_S);
        pT2_S_sum +=      pT2_S ;
        _p_thrust_pt ->fill(sqs,pT_T /MeV         );
        _p_thrust_pt2->fill(sqs,pT2_T/1e3/sqr(MeV));
        _p_sphere_pt ->fill(sqs,pT_S /MeV         );
        _p_sphere_pt2->fill(sqs,pT2_S/1e3/sqr(MeV));
      }
      _p_thrust_sum_pt ->fill(sqs,pT_T_sum /GeV);
      _p_thrust_sum_pt2->fill(sqs,pT2_T_sum/GeV);
      _p_sphere_sum_pt ->fill(sqs,pT_S_sum /GeV);
      _p_sphere_sum_pt2->fill(sqs,pT2_S_sum/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedProfilePtr<string> _p_thrust_pt, _p_thrust_pt2, _p_thrust_sum_pt, _p_thrust_sum_pt2;
    BinnedProfilePtr<string> _p_sphere_pt, _p_sphere_pt2, _p_sphere_sum_pt, _p_sphere_sum_pt2;
    string sqs;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(PLUTO_1983_I191161);

}
