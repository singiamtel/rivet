// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau -> K+ KS nu_tau
  class BABAR_2018_I1679886 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2018_I1679886);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(UnstableParticles(), "UFS");
      book(_h_KK, 1, 1, 1);
      book(_v_KK, 1, 1, 2);
      book(_c,"TMP/c");
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           unsigned int & nK0, unsigned int & nKp,
			   unsigned int & nKm, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS ) {
	  ++nKp;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::KMINUS ) {
	  ++nKm;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::K0S) {
          ++nK0;
          ++nstable;
	  ptot += p.momentum();
        }
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, nK0, nKp, nKm, ptot);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for(const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        unsigned int nstable(0),nK0(0),nKp(0),nKm(0);
      	FourMomentum p_tot(0,0,0,0);
        findDecayProducts(tau, nstable, nK0, nKp, nKm, p_tot);
        if (tau.pid() < 0) {
      	  swap(nKp,nKm);
      	}
       	if(nstable!=3) continue;
      	if(nKm==1 && nK0==1 ) {
          double q = p_tot.mass();
          _h_KK->fill(q);
          double mtau = tau.mass();
          double Cq   = q*sqr(sqr(mtau)-sqr(q))*(sqr(mtau)+2.*sqr(q));
          double fact = pow(mtau,8)/Cq;
          _v_KK->fill(q,fact);
          _c->fill();
        }
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK);
      // values from PDG 2023 for constants, see eqn 1 of paper
      double brKK = 0.5*1.486e-3; // PDG is K0 divide by 2 to get KS0
      double brE  = 0.1782;
      double Vud  = .97367;
      scale(_v_KK, 1000.*brKK/brE/sqr(Vud)/12./M_PI/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK,_v_KK;
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2018_I1679886);


}
