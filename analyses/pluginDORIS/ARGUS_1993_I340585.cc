// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D* l nu
  class ARGUS_1993_I340585 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1993_I340585);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(411); BB.addStable(-411);
      BB.addStable(421); BB.addStable(-421);
      BB.addStable(413); BB.addStable(-413);
      BB.addStable(423); BB.addStable(-423);
      BB.addStable(PID::PI0);
      declare(BB, "BB");
      // histos
      book(_h, 1,1,1);
      book(_c, "TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over B mesons
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        _c->fill();
      	for (unsigned int il=0; il < 4; ++il) {
          if (BB.modeMatches(ix,3,mode[il])) {
            int sign = il%2==0 ? -1 : 1;
            const Particle& DD = BB.decayProducts()[ix].at(sign*413)[0];
            FourMomentum qq = BB.decaying()[ix].mom() - DD.mom();
            _h->fill(qq.mass2());
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 0.5/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _c;
    const map<PdgId,unsigned int> mode[4] = {{ {-413,1}, {-11,1}, { 12,1}},
                                             { { 413,1}, { 11,1}, {-12,1}},
                                             { {-413,1}, {-13,1}, { 14,1}},
                                             { { 413,1}, { 13,1}, {-14,1}}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1993_I340585);

}
