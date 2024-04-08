// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Charged particle multiplicities and distributions
  class TASSO_1980_I143691 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1980_I143691);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Thrust(cfs), "Thrust");


      // Book histograms
      book(_mult, 1,1,1);
      unsigned int iloc(0);
      sqs = 1.0;
      if(isCompatibleWithSqrtS(13*GeV)) {
	iloc = 1;
	sqs = 13.0;
        sloc="13.0";
      }
      else if(isCompatibleWithSqrtS(17*GeV)) {
	iloc = 2;
        sqs = 17.0;
        sloc=  "17.0";
      }
      else if (isCompatibleWithSqrtS(22*GeV)) {
	iloc = 2;
        sloc= "22.0";
	sqs = 22.;
      }
      else if(isCompatibleWithSqrtS(27.6*GeV)) {
	iloc = 3;
	sqs = 27.6;
        sloc="27.6";
      }
      else if (isCompatibleWithSqrtS(30.3*GeV)) {
	iloc = 3;
	sqs = 30.3;
        sloc= "30.3";
      }
      else if (isCompatibleWithSqrtS(31.2*GeV)) {
	iloc = 3;
	sqs = 31.2;
        sloc= "31.2";
      }
      else
	MSG_ERROR("Beam energy not supported!");

      book(_h_rap,iloc+1,1,1);
      book(_h_x  ,iloc+4,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      Vector3 axis=thrust.thrustAxis();
      _mult->fill(sloc,cfs.particles().size());
      for (const Particle& p : cfs.particles()) {
	const Vector3 mom3 = p.p3();
	double pp = mom3.mod();
	double xP = 2.*pp/sqs;
	_h_x->fill(xP);
        const double mom = dot(axis, mom3);
	const double rap = 0.5 * log((p.E() + mom) /
				     (p.E() - mom));
	_h_rap->fill(fabs(rap));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_rap, 1./sumOfWeights());
      scale(_h_x  , crossSection()*sqr(sqs)/sumOfWeights()/microbarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_rap, _h_x;
    BinnedProfilePtr<string> _mult;
    double sqs;
    string sloc;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TASSO_1980_I143691);


}
