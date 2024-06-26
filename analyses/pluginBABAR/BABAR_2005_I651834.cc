// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S)/X(3782) -> J/psi pi+pi-
  class BABAR_2005_I651834 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I651834);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<int>("PID", 9030443);
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==100443 ||
						Cuts::abspid==_pid);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(443);
      declare(PSI, "PSI");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 211,1}, {-211,1}, { 443,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if(!PSI.modeMatches(ix,3,mode)) continue;
	const Particle & pip  = PSI.decayProducts()[ix].at( 211)[0];
	const Particle & pim  = PSI.decayProducts()[ix].at(-211)[0];
	double mpipi = (pip.momentum()+pim.momentum()).mass();
	_h[PSI.decaying()[ix].pid()==_pid ? 0 : 1]->fill(mpipi);
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
    int _pid;
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2005_I651834);

}
