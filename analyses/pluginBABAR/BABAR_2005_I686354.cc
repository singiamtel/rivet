// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Y(4260) -> pi+ pi_ J/psi
  class BABAR_2005_I686354 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I686354);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      _pid = getOption<int>("PID", 9030443);
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==_pid);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(443);
      declare(PSI, "PSI");
      // histograms
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 211,1}, {-211,1}, { 443,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if (!PSI.modeMatches(ix,3,mode)) continue;
	const Particle & pip  = PSI.decayProducts()[ix].at( 211)[0];
	const Particle & pim  = PSI.decayProducts()[ix].at(-211)[0];
	double mpipi = (pip.momentum()+pim.momentum()).mass();
	_h->fill(mpipi);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2005_I686354);

}
