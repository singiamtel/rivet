// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross-sections for light hadrons at 3.773 and 3.650 GeV
  class BESII_2008_I801208 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2008_I801208);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=1;ix<7;++ix) {
        stringstream ss;
        ss << "TMP/n" << ix;
        book(_nMeson[ix], 1, 1, ix);
      }
      for (const string& en : _nMeson[2].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      if(nCount[310]!=1) vetoEvent;

      if(ntotal==3) {
        if((nCount[ 211]==1 && nCount[-321]==1) ||
           (nCount[-211]==1 && nCount[ 321]==1) ) {
          if(_ecms=="3.773"s) _nMeson[1]->fill(_ecms);
        }
      }
      else if(ntotal==4&&nCount[111]==1) {
        if((nCount[ 211]==1 && nCount[-321]==1) ||
           (nCount[-211]==1 && nCount[ 321]==1) )
          _nMeson[2]->fill(_ecms);
      }
      else if(ntotal==5) {
        if((nCount[ 211]==2 && nCount[-211]==1 && nCount[-321]==1) ||
           (nCount[-211]==2 && nCount[ 211]==1 && nCount[ 321]==1) )
          _nMeson[3]->fill(_ecms);
        if(((nCount[ 211]==1 &&  nCount[-321]==1) ||
            (nCount[-211]==1 &&  nCount[ 321]==1) ) && nCount[111]==2)
          _nMeson[6]->fill(_ecms);
      }
      else if(ntotal==6&&nCount[111]==1) {
        if((nCount[ 211]==2 && nCount[-211]==1 && nCount[-321]==1) ||
           (nCount[-211]==2 && nCount[ 211]==1 && nCount[ 321]==1) )
          _nMeson[4]->fill(_ecms);
      }
      else if(ntotal==7) {
        if((nCount[ 211]==3 && nCount[-211]==2 && nCount[-321]==1) ||
           (nCount[-211]==3 && nCount[ 211]==2 && nCount[ 321]==1) ) {
          if(_ecms=="3.773"s) _nMeson[5]->fill(_ecms);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<7;++ix)
        scale(_nMeson[ix],crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nMeson[7];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2008_I801208);


}
