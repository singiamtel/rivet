// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross section for e+e- -> K+K- between 2.00 and 3.08 GeV
  class BESIII_2018_I1704558 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1704558);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nkaon, 1, 1, 1);
      for (const string& en : _nkaon.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if(abs(p.pid())!=PID::KPLUS) vetoEvent;
      }
      _nkaon->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nkaon, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nkaon;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1704558);


}
