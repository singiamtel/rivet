// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> D semileptonic
  class CLEO_1997_I443704 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1997_I443704);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(Cuts::pid==511 ||	Cuts::pid==521), "UFS");
      // histograms
      book(_h,1,1,1);
      book(_c,"TMP/nB");
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    // Calculate the recoil w using mother and daugher meson
    double recoilW(const Particle& B, int mesonID) {
      // TODO why does that not work with const?
      Particle D = select(B.children(), Cuts::pid==mesonID)[0];
      FourMomentum q = B.mom() - D.mom();
      return (B.mom()*B.mom() + D.mom()*D.mom() - q*q )/ (2. * sqrt(B.mom()*B.mom()) * sqrt(D.mom()*D.mom()) );
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double tau[2] = {1.519e-3,1.638e-3};
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
      	if (p.children().size()<=1) continue;
      	if (p.pid()==PID::BPLUS) {
      	  _c->fill();
      	  if (isSemileptonicDecay(p, {PID::D0BAR,PID::POSITRON,PID::NU_E}))  _h->fill(recoilW(p, PID::D0BAR),1./tau[1]);
      	  if (isSemileptonicDecay(p, {PID::D0BAR,PID::ANTIMUON,PID::NU_MU})) _h->fill(recoilW(p, PID::D0BAR),1./tau[1]);
      	}
      	else if (p.pid()==PID::B0) {
      	  _c->fill();
      	  if (isSemileptonicDecay(p, {PID::DMINUS,PID::POSITRON,PID::NU_E}))  _h->fill(recoilW(p, PID::DMINUS),1./tau[0]);
      	  if (isSemileptonicDecay(p, {PID::DMINUS,PID::ANTIMUON,PID::NU_MU})) _h->fill(recoilW(p, PID::DMINUS),1./tau[0]);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h,0.5/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1997_I443704);

}
