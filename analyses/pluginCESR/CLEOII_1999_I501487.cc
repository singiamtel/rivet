// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> pi- pi- pi+ pi0
  class CLEOII_1999_I501487 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1999_I501487);


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
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_2pi[ix], 1, 1, 1+ix);
        if (ix==2) continue;
        book(_h_4pi[ix], 2, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0;ix<TAU.decaying().size();++ix) {
      	int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
      	if (!(TAU.modeMatches(ix,5,mode  ) ||
              TAU.modeMatches(ix,5,modeCC))) continue;
        const Particle&  pi0 = TAU.decayProducts()[ix].at( 111     )[0];
        const Particles& pim = TAU.decayProducts()[ix].at(-211*sign);
        const Particle&  pip = TAU.decayProducts()[ix].at( 211*sign)[0];
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
        FourMomentum phad = pim[0].mom()+pim[1].mom()+pip.mom()+pi0.mom();
        const double q = phad.mass();
        if (iomega>0) {
          _h_4pi[0]->fill(q);
        }
        else {
          _h_4pi[1]->fill(q);
          _h_2pi[0]->fill((pi0.mom()+pip.mom()).mass());
          for (unsigned int iy=0; iy<2; ++iy) {
            _h_2pi[1]->fill((pim[iy].mom()+pi0.mom()).mass());
            _h_2pi[2]->fill((pim[iy].mom()+pip.mom()).mass());
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_2pi, 1.0, false);
      normalize(_h_4pi, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_2pi[3],_h_4pi[2];
    const map<PdgId,unsigned int> mode   = { { 111,1},{-211,2},{ 211,1},{ 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 111,1},{ 211,2},{-211,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1999_I501487);

}
