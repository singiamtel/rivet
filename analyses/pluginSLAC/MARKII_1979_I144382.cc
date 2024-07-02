// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> J/psi pi+pi-
  class MARKII_1979_I144382 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1979_I144382);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==100443),"UFS");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix], 1+ix, 1, 1);
    }


    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           Particles& pip, Particles& pim, Particles & onium) {
      for (const Particle& p : mother.children()) {
        int id = p.pid();
      	if (id == PID::PIMINUS) {
          pim.push_back(p);
          ++nstable;
        }
       	else if (id == PID::PIPLUS) {
       	  pip.push_back(p);
       	  ++nstable;
       	}
        else if (id==443) {
          onium.push_back(p);
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p,nstable,pip,pim,onium);
        }
        else {
          ++nstable;
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& ups : apply<UnstableParticles>(event, "UFS").particles()) {
      	unsigned int nstable(0);
      	Particles pip, pim, onium;
      	findDecayProducts(ups,nstable,pip,pim,onium);
      	// check for onium
      	if (onium.size() !=1 || nstable !=3 || pip.size()!=1 || pim.size() !=1) continue;
      	double mpipi = (pip[0].momentum()+pim[0].momentum()).mass();
      	_h[0]->fill(sqr(mpipi));
	double mpi = 2.*pip[0].mass();
	double x = (mpipi-mpi)/(ups.mass()-onium[0].mass()-mpi);
	_h[1]->fill(x);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1979_I144382);

}
