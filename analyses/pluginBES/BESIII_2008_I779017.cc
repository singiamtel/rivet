// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > light hadrons
  class BESIII_2008_I779017 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2008_I779017);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==411 ||
                                Cuts::abspid==421 ||
                                Cuts::abspid==431 ||
                                Cuts::abspid==4122||
                                Cuts::abspid==4232||
                                Cuts::abspid==4132||
                                Cuts::abspid==4332), "UFS");
      declare(FinalState(), "FS");

      // Book histograms
      book(_sigma, 1,1,1);
      _eCent = getOption<string>("ECENT", std::to_string(sqrtS()/GeV));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      if (!ufs.particles().empty()) vetoEvent;
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	vetoEvent;
      else if(nCount[-11]==1 and nCount[11]==1 &&
	      ntotal==2+nCount[22])
	vetoEvent;
      // everything else
      else
	_sigma->fill(_eCent);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma,crossSection()/sumOfWeights()/nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    string _eCent;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2008_I779017);

}
