// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambdab0 -> Lambda(1520) mu+mu-
  class LHCB_2023_I2633007 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2023_I2633007);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==5122);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 102134);
      BB.addStable(-102134);
      declare(BB, "BB");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],1,1+ix,1);
      }
      book(_c,"TMP/nLambdab");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    static const map<PdgId,unsigned int> mode   = { { 102134,1},{ 13,1}, {-13,1}};
    static const map<PdgId,unsigned int> modeCC = { {-102134,1},{ 13,1}, {-13,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        _c->fill();
      	if (! ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode)) ||
               (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,modeCC)))) continue;
        const Particle& lp = BB.decayProducts()[ix].at(-13)[0];
        const Particle& lm = BB.decayProducts()[ix].at( 13)[0];
        const double qq = (lp.mom()+lm.mom()).mass2();
        _h[0]->fill(qq);
        _h[1]->fill(qq);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1e8/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2023_I2633007);

}
