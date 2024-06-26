// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B_s -> psi(2S)/X(3782) pi+pi-
  class LHCB_2023_I2635083 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2023_I2635083);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<int>("PID", 9030443);
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==531);
      declare(ufs, "UFS");
      DecayedParticles BS0(ufs);
      BS0.addStable(100443);
      BS0.addStable(_pid);
      declare(BS0, "BS0");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    static const map<PdgId,unsigned int> mode1 = { { 211,1}, {-211,1}, { 100443,1}};
    static const map<PdgId,unsigned int> mode2 = { { 211,1}, {-211,1}, { _pid  ,1}};
      DecayedParticles BS0 = apply<DecayedParticles>(event, "BS0");
      for (unsigned int ix=0; ix<BS0.decaying().size(); ++ix) {
        unsigned int imode=0;
        if      (BS0.modeMatches(ix,3,mode1)) imode=0;
        else if (BS0.modeMatches(ix,3,mode2)) imode=1;
        else continue;
      	const Particle& pip  = BS0.decayProducts()[ix].at( 211)[0];
      	const Particle& pim  = BS0.decayProducts()[ix].at(-211)[0];
        const double mpipi = (pip.mom()+pim.mom()).mass();
      	_h[imode]->fill(mpipi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2023_I2635083);

}
