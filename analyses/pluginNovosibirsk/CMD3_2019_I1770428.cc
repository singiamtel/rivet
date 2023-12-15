// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief CMD3 K0 K0 pi+pi-
  class CMD3_2019_I1770428 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2019_I1770428);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_cK0K0pippim , "TMP/K0K0pippim");
      book( _h_pipi ,2,1,1);
      book( _h_total,2,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      Particles pip,k0;
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	if(p.abspid()==211)
	  pip.push_back(p);
	else if(p.pid()==310)
	  k0.push_back(p);
	++ntotal;
      }
      if(ntotal==4 && nCount[310]==2 && nCount[211]==1 && nCount[-211]==1) {
	_cK0K0pippim->fill();
	FourMomentum ppipi = pip[0].momentum()+pip[1].momentum();
	_h_pipi->fill(ppipi.mass()/MeV);
	for(unsigned int ix=0;ix<2;++ix)
	  _h_total->fill((ppipi+k0[ix].momentum()).mass()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi );
      normalize(_h_total);
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sigma = _cK0K0pippim->val()*fact;
      double error = _cK0K0pippim->err()*fact;
      Estimate1DPtr  mult;
      book(mult, 1, 1, 6);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/MeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _cK0K0pippim;
    Histo1DPtr _h_pipi, _h_total;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2019_I1770428);


}
