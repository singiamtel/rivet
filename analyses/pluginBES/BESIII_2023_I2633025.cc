// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta -> 3 pions
  class BESIII_2023_I2633025 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2633025);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==PID::ETA);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      ETA.addStable(PID::K0S);
      declare(ETA, "ETA");
      vector<double> bins = {-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
                             0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_d[ix],bins);
        book(_c[ix],"TMP/c_"+toString(ix));
        for(unsigned int iy=0;iy<19;++iy) {
          book(_d[ix]->bin(iy+1),3+ix,1,1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for (unsigned int ix=0; ix<ETA.decaying().size(); ++ix) {
        if (ETA.modeMatches(ix,3,mode1)) {
          const Particle& pi0 = ETA.decayProducts()[ix].at( 111)[0];
          const Particle& pip = ETA.decayProducts()[ix].at( 211)[0];
          const Particle& pim = ETA.decayProducts()[ix].at(-211)[0];
          double s1 = (pi0.mom()+pim.mom()).mass2();
          double s2 = (pi0.mom()+pip.mom()).mass2();
          double s3 = (pip.mom()+pim.mom()).mass2();
          double mOut = pi0.mass()+pip.mass()+pim.mass();
          double Q = ETA.decaying()[ix].mass()-mOut;
          double X = sqrt(3.)/2./ETA.decaying()[ix].mass()/Q*(s1-s2);
          double Y = 3.*(sqr(ETA.decaying()[ix].mass()-pi0.mass())-s3)/2./ETA.decaying()[ix].mass()/Q-1.;
          _d[0]->fill(Y,X);
          _c[0]->fill();
        }
        else if (ETA.modeMatches(ix,3,mode2)) {
          const Particles& pi0 = ETA.decayProducts()[ix].at(111);
          double s1 = (pi0[2].mom()+pi0[1].mom()).mass2();
          double s2 = (pi0[2].mom()+pi0[0].mom()).mass2();
          double s3 = (pi0[0].mom()+pi0[1].mom()).mass2();
          double mOut = pi0[2].mass()+pi0[0].mass()+pi0[1].mass();
          double Q = ETA.decaying()[ix].mass()-mOut;
          double X = sqrt(3.)/2./ETA.decaying()[ix].mass()/Q*(s1-s2);
          double Y = 3.*(sqr(ETA.decaying()[ix].mass()-pi0[2].mass())-s3)/2./ETA.decaying()[ix].mass()/Q-1.;
          _d[1]->fill(Y,X);
          _c[1]->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_d[ix],1./ *_c[ix]);
        divByGroupWidth(_d[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _d[2];
    CounterPtr _c[2];
    const map<PdgId,unsigned int> mode1 = { {211,1}, {-211,1}, {111,1} };
    const map<PdgId,unsigned int> mode2 = { {111,3} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2633025);

}
