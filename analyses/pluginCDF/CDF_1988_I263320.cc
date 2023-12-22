// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerCDFRun0Run1.hh"

namespace Rivet {


  /// @brief CDF track \f$ p_\perp \f$ distributions at 630 and 1800 GeV
  class CDF_1988_I263320 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1988_I263320);


    /// @name Analysis methods
    /// @{

    /// Book histograms and set up projections
    void init() {
      // Set up projections
      declare(TriggerCDFRun0Run1(), "Trigger");
      const ChargedFinalState cfs((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  0.4*GeV));
      declare(cfs, "CFS");

      // Book histo
      if (isCompatibleWithSqrtS(1800*GeV)) {
        book(_hist_pt, 1, 1, 1);
        _axis = YODA::Axis<double>{0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0125, 1.1, 1.2, 1.3,
                                   1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.025,
                                   3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.05, 5.4, 5.8, 6.2, 6.6, 7.15, 8.25, 9.75};
      }
      else if (isCompatibleWithSqrtS(630*GeV)) {
        book(_hist_pt, 2, 1, 1);
        _axis = YODA::Axis<double>{0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0125, 1.1, 1.2, 1.3,
                                   1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.025, 2.2, 2.4, 2.6, 2.8, 3.2, 3.8};
      }

      book(_sumWTrig, "sumWTrig");

    }


    /// Do the analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _hist_pt->xEdges();
      // Trigger
      const bool trigger = apply<TriggerCDFRun0Run1>(event, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;
      _sumWTrig->fill();

      const FinalState& trackfs = apply<ChargedFinalState>(event, "CFS");
      for (Particle p : trackfs.particles()) {
        const double pt = p.pT()/GeV;
        // Effective weight for d3sig/dp3 = weight / ( Delta eta * 2pi * pt ), with Delta(eta) = 2
        const double eff_weight = 1.0/(2*2*TWOPI*pt);
        _hist_pt->fill(map2string(pt), eff_weight);
      }
    }


   string map2string(const double value) const {
     const size_t idx = _axis.index(value);
     if (idx && idx < _edges.size() + 1)  return _edges[idx-1];
     return "OTHER";
   }

    /// Scale histos
    void finalize() {
      scale(_hist_pt, crossSectionPerEvent()/millibarn);
      for (auto& b : _hist_pt->bins()) {
        b.scaleW(1.0/_axis.width(b.index()));
      }
    }

    /// @}


  private:

    /// Counter
    CounterPtr _sumWTrig;

    /// Histo
    BinnedHistoPtr<string> _hist_pt;
    YODA::Axis<double> _axis;
    vector<string> _edges;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_1988_I263320, CDF_1988_S1865951);

}
