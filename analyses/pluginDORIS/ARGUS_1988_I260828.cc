// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> K+ K- pi+ pi-
  class ARGUS_1988_I260828 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1988_I260828);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      book(_h, 1, 1, 1);

      _edge = "OTHER"s;
      for (const string& edge : _h.binning().edges<0>()) {
        if (isCompatibleWithSqrtS(std::stod(edge)*GeV)) {
          _edge = edge;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      if (ntotal==5 && nCount[211]==1 && nCount[-211]==1 &&
          nCount[321]==1 && nCount[-321]==1 && nCount[111]==1) _h->fill(_edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection()/nanobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h;
    string _edge;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1988_I260828);

}
