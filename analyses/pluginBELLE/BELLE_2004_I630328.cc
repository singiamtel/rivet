// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> p pbar pi+/K+/K*+ and B0 -> p pbar KS0
  class BELLE_2004_I630328 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I630328);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 310);
      BB.addStable( 323);
      BB.addStable(-323);
      declare(BB, "BB");
      // histograms
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_c[ix],"TMP/nB_"+toString(ix+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int imode=0;
        if (BB.decaying()[ix].abspid()==511) _c[0]->fill();
        else                                 _c[1]->fill();
        if ((BB.decaying()[ix].pid()   == 521 && BB.modeMatches(ix,3,mode1)) ||
            (BB.decaying()[ix].pid()   ==-521 && BB.modeMatches(ix,3,mode1CC))) imode = 0;
        else if ((BB.decaying()[ix].pid()   == 521 && BB.modeMatches(ix,3,mode2)) ||
                 (BB.decaying()[ix].pid()   ==-521 && BB.modeMatches(ix,3,mode2CC))) imode = 1;
        else if ((BB.decaying()[ix].pid()   == 521 && BB.modeMatches(ix,3,mode3)) ||
                 (BB.decaying()[ix].pid()   ==-521 && BB.modeMatches(ix,3,mode3CC))) imode = 2;
        else if (BB.decaying()[ix].abspid()== 511 && BB.modeMatches(ix,3,mode4))  imode = 3;
        else {
          continue;
        }
        const Particle& pp   = BB.decayProducts()[ix].at( 2212)[0];
        const Particle& pbar = BB.decayProducts()[ix].at(-2212)[0];
        const double mass = (pp.mom()+pbar.mom()).mass();
       	_h[imode]->fill(mass);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<4; ++ix) {
        if (ix==2) scale(_h[ix], 1e6/ *_c[0]);
        else       scale(_h[ix], 1e6/ *_c[1]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    CounterPtr _c[2];
    const map<PdgId,unsigned int> mode1   = { { 2212,1}, {-2212,1}, { 321,1} };
    const map<PdgId,unsigned int> mode1CC = { { 2212,1}, {-2212,1}, {-321,1} };
    const map<PdgId,unsigned int> mode2   = { { 2212,1}, {-2212,1}, { 211,1} };
    const map<PdgId,unsigned int> mode2CC = { { 2212,1}, {-2212,1}, {-211,1} };
    const map<PdgId,unsigned int> mode3   = { { 2212,1}, {-2212,1}, { 323,1} };
    const map<PdgId,unsigned int> mode3CC = { { 2212,1}, {-2212,1}, {-323,1} };
    const map<PdgId,unsigned int> mode4   = { { 2212,1}, {-2212,1}, { 310,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2004_I630328);

}
