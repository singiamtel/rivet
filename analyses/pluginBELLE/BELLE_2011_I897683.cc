// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> phi K gamma
  class BELLE_2011_I897683 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2011_I897683);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PHI);
      BB.addStable(PID::K0S);
      declare(BB, "BB");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
      	int sign = 1, imode=0, iK=321;
        if (BB.decaying()[ix].abspid()==521) {
          if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
            sign=1;
          }
          else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
            sign=-1;
          }
          else {
            continue;
          }
        }
        else if (BB.modeMatches(ix,3,mode2)) {
          imode=1;
          iK=310;
        }
        else {
          continue;
        }
        const Particle& K   = BB.decayProducts()[ix].at( sign*iK)[0];
        const Particle& phi = BB.decayProducts()[ix].at(     333)[0];
        _h[imode]->fill((K.mom()+phi.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1., false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    const map<PdgId,unsigned int> mode1   = { { 333,1},{ 321,1}, {22,1}};
    const map<PdgId,unsigned int> mode1CC = { { 333,1},{-321,1}, {22,1}};
    const map<PdgId,unsigned int> mode2   = { { 333,1},{ 310,1}, {22,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2011_I897683);

}
