// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief b-baryon polarization
  class ALEPH_1996_I402895 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_1996_I402895);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_El, "El",   45,0.,45.0);
      book(_h_Ev, "Ev",   45,0.,45.0);

    }

    void findDecayProducts(Particle p, Particles & lep, Particles & nu) {
      for(const Particle & child : p.children()) {
        if(PID::isHadron(child.pid())) continue;
        if(child.abspid()==11 or child.abspid()==13)
          lep.push_back(child);
        else if(child.abspid()==12 or child.abspid()==14)
          nu.push_back(child);
        else if(child.abspid()!=15)
          findDecayProducts(child,lep,nu);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over weakly decaying b-baryons
      for (const Particle& p : ufs.particles(Cuts::abspid==5122 || Cuts::abspid==5132 ||
					     Cuts::abspid==5232 || Cuts::abspid==5332)) {
        Particles lep,nu;
        findDecayProducts(p,lep,nu);
        if(lep.size()!=1 || nu.size()!=1) continue;
        _h_El   ->fill(lep[0].momentum().t());
        _h_Ev   ->fill(nu [0].momentum().t());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_El   );
      normalize(_h_Ev   );
      if(_h_El->effNumEntries()!=0. and _h_Ev->effNumEntries()!=0.) {
        const double Ev  = _h_Ev->xMean();
        const double El  = _h_El->xMean();
        const double dEv = _h_Ev->xStdErr();
        const double dEl = _h_El->xStdErr();
        const double ratio = El/Ev;
        const double dr    = (Ev*dEl-El*dEv)/sqr(Ev);
        const double rho = 0.091;
        const double P = (7. + rho*(30. - 40.*ratio) + 4.*(2.-3.*ratio)*ratio)/sqr(1.+2.*ratio);
        const double dP = (20.*(-1. + 4.*rho*(-2. + ratio) - 2.*ratio))/pow(1. + 2.*ratio,3)*dr;
       	BinnedEstimatePtr<string> h_pol;
        book(h_pol,1,1,1);
       	h_pol->bin(1).set(P, dP);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_El, _h_Ev;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALEPH_1996_I402895);


}
