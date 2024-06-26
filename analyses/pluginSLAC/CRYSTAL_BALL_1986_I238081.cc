// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > hadrons
  class CRYSTAL_BALL_1986_I238081 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1986_I238081);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons",refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons,   "/TMP/sigma_muons"  ,refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_D_star,  "/TMP/sigma_D_star" ,refData(2,1,1));
      
      for (const string& ecms : _c_hadrons.binning().edges<0>()) {
        if (isCompatibleWithSqrtS(std::stod(ecms)*GeV)) {
          _ecms = ecms;
          break;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // mu+mu- + photons
      if(!_ecms.empty()) {
        if(nCount[-13]==1 and nCount[13]==1 &&
           ntotal==2+nCount[22])
          _c_muons->fill(_ecms);
        // everything else
        else
          _c_hadrons->fill(_ecms);
      }
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      bool found = false;
      for (const Particle & p : ufs.particles()) {
	if(abs(p.pid())!=413 && abs(p.pid())!=423) continue;
	bool fs = true;
	for (const Particle & child : p.children()) {
	  if(child.pid()==p.pid()) {
	    fs = false;
	    break;
	  }
	}
	if(fs) {
	  found = true;
	  break;
	}
      }
      if(found)
	_c_D_star->fill(sqrtS());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // R
      BinnedEstimatePtr<string> tmp;
      book(tmp,1,1,1);
      divide(_c_hadrons,_c_muons,tmp);
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale(_c_D_star,fact);
      Estimate1DPtr tmp2;
      book(tmp2,2,1,1);
      barchart(_c_D_star,tmp2);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    Histo1DPtr _c_D_star;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1986_I238081);


}
