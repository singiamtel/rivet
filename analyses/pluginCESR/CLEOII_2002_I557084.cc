// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> D*0 2pi+2pi-
  class CLEOII_2002_I557084 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2002_I557084);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::PI0);
      B0.addStable( 413);
      B0.addStable(-413);
      B0.addStable( 423);
      B0.addStable(-423);
      B0.addStable( 411);
      B0.addStable(-411);
      B0.addStable( 421);
      B0.addStable(-421);
      B0.addStable( 431);
      B0.addStable(-431);
      declare(B0, "B0");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
       	book(_h_Dpi[ix],1,1,1+ix);
      	book(_h_4pi[ix],3,1,1+ix);
      }
      book(_h_3pi,2,1,1);
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over particles
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      int imode = -1;
      for (unsigned int ix=0; ix<B0.decaying().size(); ++ix) {
        int sign = B0.decaying()[ix].pid()/B0.decaying()[ix].abspid();
        _c->fill();
        if ( (sign== 1 && B0.modeMatches(ix,5,mode1) ) ||
             (sign==-1 && B0.modeMatches(ix,5,mode1CC) ) ) {
          imode=0;
        }
        else if ( (sign== 1 && B0.modeMatches(ix,5,mode2) ) ||
                  (sign==-1 && B0.modeMatches(ix,5,mode2CC) ) ) {
          imode=1;
        }
        else {
          continue;
        }
        const Particles& pip = B0.decayProducts()[ix].at( sign*211);
        const Particles& pim = B0.decayProducts()[ix].at(-sign*211);
        FourMomentum pTotal;
        for (const Particle& p : pip) pTotal+=p.mom();
        for (const Particle& p : pim) pTotal+=p.mom();
        if (imode==0) {
          const Particle& pi0 = B0.decayProducts()[ix].at(111)[0];
          pTotal+= pi0.mom();
          _h_4pi[1]->fill(pTotal.mass2());
          continue;
        }
        _h_4pi[0]->fill(pTotal.mass2());
        // D* pi mases
        const Particle& Dstar = B0.decayProducts()[ix].at(-sign*423)[0];
        for (const Particle& p : pim) {
          FourMomentum pDpi = Dstar.mom()+p.mom();
          const double mDpi = pDpi.mass();
          _h_Dpi[0]->fill(mDpi);
          if (mDpi>2.3 and mDpi<2.6) {
            FourMomentum p3pi = B0.decaying()[ix].mom()-pDpi;
            _h_3pi->fill(p3pi.mass());
          }
        }
        for (const Particle& p : pip) {
          _h_Dpi[1]->fill((Dstar.mom()+p.mom()).mass());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // first hist is differential BR/ BR for leptons
      const double br = .1033;
      normalize(_h_Dpi, 1.0, false);
      scale(_h_4pi, 1.0/br/ *_c);
      normalize(_h_3pi,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Dpi[2], _h_3pi, _h_4pi[2];
    CounterPtr _c;
    const map<PdgId,unsigned int> mode1   = { { -413,1}, { 211,2}, {-211,1}, {111,1}};
    const map<PdgId,unsigned int> mode1CC = { {  413,1}, {-211,2}, { 211,1}, {111,1}};
    const map<PdgId,unsigned int> mode2   = { { -423,1}, { 211,2}, {-211,2}};
    const map<PdgId,unsigned int> mode2CC = { {  423,1}, { 211,2}, {-211,2}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2002_I557084);

}
