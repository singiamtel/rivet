// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- > pi+ pi- 2pi0 (including omega)
  class BABAR_2017_I1621593 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2017_I1621593);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix)
        book(_mult[ix], "TMP/mult_"+toString(ix), refData(1+ix, 1, 1));
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
      if(ntotal!=4) vetoEvent;
      if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==2) {
        _mult[0]->fill(sqrtS()/GeV);
        const FinalState& ufs = apply<FinalState>(event, "UFS");
        if (!ufs.particles(Cuts::pid==223).empty()) {
          _mult[1]->fill(sqrtS()/GeV);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(size_t ix=0;ix<2;++ix) {
        scale(_mult[ix],fact);
        Estimate1DPtr tmp;
        book(tmp,1+ix,1,1);
        barchart(_mult[ix],tmp);
      }
    }

    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _mult[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2017_I1621593);


}
