// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> 3(pi+pi-)
  class CMD3_2013_I1217420 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2013_I1217420);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms

      book(_c_all, "/TMP/all");

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the final-state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==6 && nCount[211]==3 && nCount[-211]==3) {
	_c_all->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact = crossSection()/nanobarn/sumOfWeights();
      double sigma = _c_all->val()*fact;
      double error = _c_all->err()*fact;
      for (unsigned int ihist=1;ihist<4;++ihist) {
        Estimate1DPtr  mult;
        book(mult, 1, 1, ihist);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/MeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_all;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2013_I1217420);


}
