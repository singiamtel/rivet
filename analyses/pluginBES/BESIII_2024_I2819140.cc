// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> K- ell+ nu_ell and D+ -> Kbar0 ell+ nu_ell
  class BESIII_2024_I2819140 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2819140);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==411 || Cuts::abspid==421 ), "UFS");

      for(unsigned int ix=0;ix<4;++ix) {
        book(_h[ix], 1+ix , 1 ,1);
        if(ix>1) continue;
        book(_c[ix], "TMP/c_"+toString(ix+1));
      }
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - select(B.children(), Cuts::pid==mesonID)[0];
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over D mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        if(p.pid()==421) {
          _c[0]->fill();
          if(isSemileptonicDecay(p, {PID::KMINUS, PID::EPLUS, PID::NU_E})) {
            _h[0]->fill(q2(p, PID::KMINUS));
          }
          else if(isSemileptonicDecay(p, {PID::KMINUS, PID::ANTIMUON, PID::NU_MU})) {
            _h[1]->fill(q2(p, PID::KMINUS));
          }
        }
        else if(p.pid()==-421) {
          _c[0]->fill();
          if(isSemileptonicDecay(p, {PID::KPLUS, PID::EMINUS, PID::NU_EBAR})) {
            _h[0]->fill(q2(p, PID::KPLUS));
          }
          else if(isSemileptonicDecay(p, {PID::KPLUS, PID::MUON, PID::NU_MUBAR})) {
            _h[1]->fill(q2(p, PID::KPLUS));
          }
        }
        else if(p.pid()==411) {
          _c[1]->fill();
          if(isSemileptonicDecay(p, {-311, PID::EPLUS, PID::NU_E})) {
            _h[2]->fill(q2(p, -311));
          }
          else if(isSemileptonicDecay(p, {-311, PID::ANTIMUON, PID::NU_MU})) {
            _h[3]->fill(q2(p, -311));
          }
          else if(isSemileptonicDecay(p, { 130, PID::EPLUS, PID::NU_E})) {
            _h[2]->fill(q2(p, 130));
          }
          else if(isSemileptonicDecay(p, { 130, PID::ANTIMUON, PID::NU_MU})) {
            _h[3]->fill(q2(p, 130));
          }
          else if(isSemileptonicDecay(p, { 310, PID::EPLUS, PID::NU_E})) {
            _h[2]->fill(q2(p, 310));
          }
          else if(isSemileptonicDecay(p, { 310, PID::ANTIMUON, PID::NU_MU})) {
            _h[3]->fill(q2(p, 310));
          }
        }
        else if(p.pid()==-411) {
          _c[1]->fill();
          if(isSemileptonicDecay(p, {311, PID::EMINUS, PID::NU_EBAR})) {
            _h[2]->fill(q2(p, 311));
          }
          else if(isSemileptonicDecay(p, {311, PID::MUON, PID::NU_MUBAR})) {
            _h[3]->fill(q2(p, 311));
          }
          else if(isSemileptonicDecay(p, { 130, PID::EMINUS, PID::NU_EBAR})) {
            _h[2]->fill(q2(p, 130));
          }
          else if(isSemileptonicDecay(p, { 130, PID::MUON, PID::NU_MUBAR})) {
            _h[3]->fill(q2(p, 130));
          }
          else if(isSemileptonicDecay(p, { 310, PID::EMINUS, PID::NU_EBAR})) {
            _h[2]->fill(q2(p, 310));
          }
          else if(isSemileptonicDecay(p, { 310, PID::MUON, PID::NU_MUBAR})) {
            _h[3]->fill(q2(p, 310));
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // lifetimes from PDG 2022 [ns] (as in paper)
      double tau[2] = {410.3e-6, 1033e-6};
      for(unsigned int ix=0;ix<4;++ix)
        scale(_h[ix],1./tau[ix/2]/ *_c[ix/2]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2819140);

}
