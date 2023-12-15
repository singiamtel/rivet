// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief CDF Z boson rapidity measurement
  class CDF_2009_I856131 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2009_I856131);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      // this seems to have been corrected completely for all selection cuts,
      // i.e. eta cuts and pT cuts on leptons.
      DileptonFinder zfinder(91.2*GeV, 0.2, Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
      declare(zfinder, "DileptonFinder");


      /// Book histograms here
      book(_h_xs ,1, 1, 1);
      book(_h_yZ ,2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const DileptonFinder& zfinder = apply<DileptonFinder>(event, "DileptonFinder");
      if (zfinder.bosons().size() == 1) {
        _h_yZ->fill(fabs(zfinder.bosons()[0].rapidity()));
        _h_xs->fill(1960);
      } else {
        MSG_DEBUG("no unique lepton pair found.");
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_xs, crossSection()/sumOfWeights());
      // Data seems to have been normalized for the avg of the two sides
      // (+ve & -ve rapidity) rather than the sum, hence the 0.5:
      scale(_h_yZ, 0.5*crossSection()/sumOfWeights());
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_yZ;
    Histo1DPtr _h_xs;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(CDF_2009_I856131);

}
