// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> pi, K ell nu_ell
  class BELLE_2006_I715430 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I715430);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==PID::D0), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],"TMP/h_"+toString(ix+1), refData(1+ix, 1, 1));
      }
      book(_c,"TMP/nD");
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) const {
      FourMomentum q = B.mom() - select(B.children(), Cuts::abspid==abs(mesonID))[0];
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) const {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double M0 = 1.86484, M02 = sqr(M0);
      const double m2[2] = {sqr(0.493677),sqr(0.13957039)};
      // Loop over D mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        _c->fill();
        if (isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E}) ||
            isSemileptonicDecay(p, {PID::PIPLUS , PID::ELECTRON, PID::NU_EBAR}) ) {
          double qq = q2(p, PID::PIMINUS);
          double pcm = sqrt(0.25/M02*(sqr(M02)+sqr(qq)+sqr(m2[1])-2.*qq*m2[1]-2.*M02*qq-2.*m2[1]*M02));
          _h[1]->fill(qq,1./pow(pcm,3));
        }
        else if(isSemileptonicDecay(p, {PID::KMINUS, PID::POSITRON, PID::NU_E}) ||
                isSemileptonicDecay(p, {PID::KPLUS , PID::ELECTRON, PID::NU_EBAR})) {
          double qq = q2(p, PID::KMINUS);
          double pcm = sqrt(0.25/M02*(sqr(M02)+sqr(qq)+sqr(m2[0])-2.*qq*m2[0]-2.*M02*qq-2.*m2[0]*M02));
          _h[0]->fill(qq,1./pow(pcm,3));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double GF      = 1.1663788e-5;
      const double CKM[2]  = {0.975,0.221};
      const double gamma   = 6.582119569e-25/4.103e-13;
      const double pre    = 24.*pow(M_PI,3)/sqr(GF)*gamma;
      for (unsigned int ix=0; ix<2; ++ix) {
        double fact = pre/sqr(CKM[ix]);
        Estimate1DPtr tmp;
        book(tmp,1+ix,1,1);
        scale(_h[ix],fact/ *_c);
        for (unsigned int iy=0; iy<_h[ix]->numBins(); ++iy) {
          double val = _h[ix]->bin(iy+1).sumW()/_h[ix]->bin(iy+1).xWidth();
          double err = _h[ix]->bin(iy+1).errW()/_h[ix]->bin(iy+1).xWidth();
          val = sqrt(val);
          err /= 2.*val ;
          tmp->bin(iy+1).set(val,err);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2006_I715430);

}
