// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief An example of using HDF5 analysis aux data
  class EXAMPLE_AUXDATA : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_AUXDATA);


    /// @name Analysis methods
    /// @{

    /// Initialise projections and output before the run
    void init() {

      // Projections
      FinalState fs(Cuts::abseta < 4.9);
      FastJets jetfs(fs, JetAlg::ANTIKT, 0.4);
      declare(jetfs, "Jets");
      
      // Get scale-factors from aux data
      auxData("jetpts", _jetpts);
      auxData("jetscales", _jetscales);

      // Histogram
      book(_h_jetmass, "jetmass", logspace(10, 20., 100.));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve objects
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV);

      // Fill each jet mass into a histogram with a pT-indexed reco scale-factor
      for (const Jet& j : jets) {
	double sf = _jetscales[ binIndex(j.pT()/GeV, _jetpts, true) ];
	_h_jetmass->fill(j.mass()/GeV * sf);
      }

    }

    
    /// Finalize (nothing to do since ofstream closes automatically)
    // void finalize() {}

    /// @}


    /// Scalefactors
    vector<double> _jetpts, _jetscales;

    // Histogram
    Histo1DPtr _h_jetmass;
    
  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_AUXDATA);

}
