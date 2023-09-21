// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross-sections for light hadrons at 3.773, 3.650 and 3.6648 GeV
  class BESII_2008_I801210 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2008_I801210);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=1;ix<6;++ix) {
        stringstream ss;
        ss << "TMP/n" << ix;
        book(_nMeson[ix], ss.str());
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
      if(nCount[111]!=2) vetoEvent;

      if(ntotal==4) {
        if(nCount[211]==1 && nCount[-211]==1)
          _nMeson[1]->fill();
        else if(nCount[321]==1 && nCount[-321]==1)
          _nMeson[2]->fill();
      }
      else if(ntotal==6) {
        if(nCount[211]==2 && nCount[-211]==2)
          _nMeson[3]->fill();
        else if(nCount[321]==1 && nCount[-321]==1 &&
                nCount[211]==1 && nCount[-211]==1)
          _nMeson[4]->fill();
      }
      else if(ntotal==8) {
        if(nCount[211]==3 && nCount[-211]==3)
          _nMeson[5]->fill();
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<6;++ix) {
        double sigma = _nMeson[ix]->val();
        double error = _nMeson[ix]->err();
        sigma *= crossSection()/ sumOfWeights() /picobarn;
        error *= crossSection()/ sumOfWeights() /picobarn;
        Estimate1DPtr  mult;
        book(mult, 1, 1, ix);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[6];
    /// @}


  };



  RIVET_DECLARE_PLUGIN(BESII_2008_I801210);

}
