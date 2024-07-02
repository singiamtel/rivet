// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief omega -> pi0 e+e-
  class CMD2_2005_I676548 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_2005_I676548);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==223), "UFS");
      // Book histograms
      book(_h_omega, 2, 1, 1);
      book(_weight_omega,"TMP/weight_omega");
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & npi, 
                           unsigned int & nep, unsigned int & nem, unsigned int & ngamma,
			   FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::EMINUS ) {
          ++nem;
          ++nstable;
          ptot += p.momentum();
        }
        else if (id == PID::EPLUS) {
          ++nep;
          ++nstable;
          ptot += p.momentum();
        }
        else if (id == PID::PI0) {
          ++npi;
          ++nstable;
        }
        else if (id == PID::GAMMA && p.children().empty() ) {
          ++ngamma;
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, npi,nep,nem,ngamma,ptot);
        }
        else {
          ++nstable;
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double me     = 0.5109989461*MeV;
      const double momega = 782.65*MeV;
      const double mpi    = 134.9770*MeV;
      // Loop over eta and omega mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
       	unsigned nstable(0),npi(0),nep(0),nem(0),ngamma(0);
       	FourMomentum ptot;
      	findDecayProducts(p,nstable,npi,nep,nem,ngamma,ptot);
        if (nstable==3 && nem==1 && nep==1 && npi==1) {
          const double q = ptot.mass();
          const double beta = sqrt(1.-4*sqr(me/q));
          const double p = sqrt(sqr(1.+sqr(q)/(sqr(momega)-sqr(mpi)))-4.*sqr(momega*q/(sqr(momega)-sqr(mpi))));
          const double fact = beta*GeV/q*(1.+2.*sqr(me/q))*pow(p,3);
          _h_omega->fill(q/GeV,1./fact);
        }
        else if(nstable==2 && ngamma ==1 && npi==1) {
          _weight_omega->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double alpha= 7.2973525664e-3;
      scale(_h_omega, 1.5 *M_PI/alpha/ *_weight_omega);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_omega;
    CounterPtr _weight_omega;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD2_2005_I676548);

}
