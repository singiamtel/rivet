// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> 4 pions
  class BESIII_2024_I2730532 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2730532);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      for(unsigned int ix=0;ix<7;++ix) {
        book(_h_0[ix],2,1,1+ix);
        if(ix>=3) continue;
        book(_h_c[ix],1,1,1+ix);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for (unsigned int ix=0;ix<D0.decaying().size(); ++ix) {
      int sign = D0.decaying()[ix].pid()/421;
      if (D0.modeMatches(ix,4,mode1)) {
        const Particles & pip= D0.decayProducts()[ix].at( sign*211);
        const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
        bool Kveto=false;
        set<double> mpm;
        for (unsigned int ix=0;ix<2;++ix) {
          if((pim[0].mom()+pim[1].mom()+pip[ix].mom()).mass()<0.51) continue;
          for(unsigned int iy=0;iy<2;++iy) {
            double m = (pip[ix].mom()+pim[iy].mom()).mass();
            mpm.insert(m);
            if(abs(m-.4976)<0.03) Kveto=true;
          }
        }
        if (Kveto) continue;
        for (const double& m : mpm) _h_c[0]->fill(m);
          for (unsigned int ix=0;ix<2;++ix) {
            _h_c[1]->fill((D0.decaying()[ix].mom()-pim[ix].mom()).mass());
            _h_c[2]->fill((D0.decaying()[ix].mom()-pip[ix].mom()).mass());
          }
        }
        else if (D0.modeMatches(ix,4,mode2)) {
          const Particles& pi0= D0.decayProducts()[ix].at(      111);
          const Particle & pip= D0.decayProducts()[ix].at( sign*211)[0];
          const Particle & pim= D0.decayProducts()[ix].at(-sign*211)[0];
          double m00 = (pi0[0].mom()+pi0[1].mom()).mass();
          if (m00>0.4376&&m00<5276) continue;
          double mpm =  (pip.mom()+pim.mom()).mass();
          if (abs(mpm-.4976)<0.03) continue;
          double mp0[2],mm0[2],mpm0[2];
          bool veto=false;
          for (unsigned int ix=0;ix<2;++ix) {
            mp0[ix] = (pi0[ix].mom()+pip.mom()).mass();
            mm0[ix] = (pi0[ix].mom()+pim.mom()).mass();
            if (mm0[ix]>0.4677&&mm0[ix]<.5067) veto=true;
            mpm0[ix] = (pi0[ix].mom()+pim.mom()+pip.mom()).mass();
            if (mpm0[ix]<.57) veto=true;
          }
          if (veto) continue;
          for (unsigned int ix=0;ix<2;++ix) {
            _h_0[0]->fill(mpm0[ix]);
            _h_0[5]->fill(mp0[ix]);
            _h_0[6]->fill(mm0[ix]);
          }
          _h_0[1]->fill((pi0[0].mom()+pi0[1].mom()+pip.mom()).mass());
          _h_0[2]->fill((pi0[0].mom()+pi0[1].mom()+pim.mom()).mass());
          _h_0[3]->fill(mpm);
          _h_0[4]->fill(m00);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_c);
      normalize(_h_0);
    }
    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _h_c[3],_h_0[7];
    const map<PdgId,unsigned int> mode1 = { { 211,2}, { -211,2}};
    const map<PdgId,unsigned int> mode2 = { { 111,2}, { 211,1}, { -211,1}};
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2730532);

}
