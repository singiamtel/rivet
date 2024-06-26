// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @brief Direct photon spectrum in upslion decay
  class CLEOIII_2006_I701217 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOIII_2006_I701217);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==553 ||
				Cuts::pid==100553 ||
				Cuts::pid==200553), "UFS");
      // histos
      for(unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix], 1+ix, 1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles()) {
      	unsigned int nhadron(0);
      	Particles photons;
      	for (const Particle & child : ups.children()) {
      	  if (PID::isHadron(child.pid())) {
      	    ++nhadron;
          }
      	  else if (child.pid()==PID::PHOTON) {
      	    photons.push_back(child);
          }
      	}
        if (nhadron==0 || photons.empty()) continue;
        LorentzTransform boost;
        if (ups.p3().mod() > 1*MeV) {
          boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
        }
        unsigned int iups = ups.pid()/100000;
        for (const Particle & gamma : photons) {
          double E = boost.transform(gamma.momentum()).E();
          double sigma = 0.01*E*(0.6*pow(E,-0.7309)+1.14-0.01*E);
	        E = randnorm(E,sigma);
          _h[iups]->fill(2.*E/ups.mass());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOIII_2006_I701217);

}
