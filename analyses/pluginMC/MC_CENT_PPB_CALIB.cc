// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Analyses/MC_CENT_PPB_Projections.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"

namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  class MC_CENT_PPB_CALIB : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(MC_CENT_PPB_CALIB);


    /// Book histograms and initialise projections before the run
    void init() {

      // One projection for the actual observable, and one for the
      // generated impact parameter.
      declare(MC_SumETFwdPbCentrality(), "Centrality");
      declare(ImpactParameterProjection(), "IMP");
      declare(MC_pPbMinBiasTrigger(), "Trigger");

      // The calibrationhistogram:
      book(_calib, "SumETPb", 100, 0.0, 200.0);

      // If histogram was pre-loaded, the calibration is done.
      _done = ( _calib->numEntries() > 0 );

      // The alternative histogram based on impact parameter. Note that
      // it MUST be named the same as the histogram for the experimental
      // observable with an added _IMP suffix for the Pecentile<>
      // binning to work properly.
      book(_impcalib, "SumETPb_IMP", 400, 0.0, 20.0);


    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if ( _done ) return;

      // The alternative centrality based on generated impact
      // parameter, assumes that the generator does not describe the
      // full final state, and should therefore be filled even if the
      // event is not triggered.
      _impcalib->fill(apply<SingleValueProjection>(event, "IMP")());

      if ( !apply<TriggerProjection>(event, "Trigger")() ) vetoEvent;

      _calib->fill(apply<SingleValueProjection>(event, "Centrality")());

    }

    /// Finalize
    void finalize() {

      _calib->normalize();
      _impcalib->normalize();

    }

  private:

    /// The calibration histograms.
    Histo1DPtr _calib;
    Histo1DPtr _impcalib;

    /// Safeguard from adding to a pre-loaded histogram.
    bool _done;

  };


  RIVET_DECLARE_PLUGIN(MC_CENT_PPB_CALIB);

}
