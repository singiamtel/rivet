// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of R
  class KEDR_2019_I1673357 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KEDR_2019_I1673357);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==443 or Cuts::pid==100443
                                or Cuts::pid==30443), "UFS");

      // Book histograms
      book(_c_hadrons[0], "/TMP/sigma_hadrons_1",refData(1,1,1));
      book(_c_hadrons[1], "/TMP/sigma_hadrons_2",refData(1,1,2));
      book(_c_muons     , "/TMP/sigma_muons"    ,refData(1,1,1));
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
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill(sqrtS()/MeV);
      // everything else
      else {
        const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
        if(ufs.particles().empty()) _c_hadrons[0]->fill(sqrtS()/MeV);
        _c_hadrons[1]->fill(sqrtS()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        Estimate1DPtr mult;
        book(mult, 1, 1, 1+ix);
        divide(_c_hadrons[ix],_c_muons,mult);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _c_hadrons[2], _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KEDR_2019_I1673357);

}
