// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> pi- pi- pi+ pi0
  class CLEOII_1995_I400623 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1995_I400623);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      TAU.addStable(221);
      declare(TAU, "TAU");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_spect[ix], 1, 1, 1+ix);
      }
      book(_c,"TMP/nTau");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0; ix<TAU.decaying().size(); ++ix) {
        _c->fill();
      	int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
      	if (!(TAU.modeMatches(ix,5,mode) || TAU.modeMatches(ix,5,modeCC))) continue;
        const Particle&  pi0 = TAU.decayProducts()[ix].at( 111     )[0];
        const Particles& pim = TAU.decayProducts()[ix].at(-211*sign);
        const Particle&  pip = TAU.decayProducts()[ix].at( 211*sign)[0];
        FourMomentum phad = pim[0].mom()+pim[1].mom()+pip.mom()+pi0.mom();
        const double q = phad.mass();
        const double mtau = TAU.decaying()[ix].mass();
        _h_spect[0]->fill(phad.mass(),sqr(mtau)*mtau/q/sqr(sqr(mtau)-sqr(q))/(sqr(mtau)+2.*sqr(q)));
        Particle omega;
        int iomega=-1;
        for (unsigned int ix=0; ix<2; ++ix) {
          Particle parent = pim[ix];
          while (!parent.parents().empty()) {
            parent = parent.parents()[0];
            if (parent.pid()==223) break;
          }
          if (parent.pid()==223) {
            omega=parent;
            iomega=ix;
          }
        }
        if (iomega<0) continue;
        _h_spect[1]->fill(phad.mass(),sqr(mtau)*mtau/q/sqr(sqr(mtau)-sqr(q))/(sqr(mtau)+2.*sqr(q)));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double Gamma = 6.582119514e-25/290.3e-15;
      const double GF    = 1.1663787e-5;
      const double cosC  = 0.97373;
      scale(_h_spect, Gamma*16.*sqr(M_PI)/sqr(GF*cosC)/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect[2];
    CounterPtr _c;
    const map<PdgId,unsigned int> mode   = { { 111,1},{-211,2},{ 211,1},{ 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 111,1},{ 211,2},{-211,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1995_I400623);

}
