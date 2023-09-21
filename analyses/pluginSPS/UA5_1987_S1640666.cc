// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TriggerUA5.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// UA5 charged multiplicity measurements at 546 GeV
  class UA5_1987_S1640666 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(UA5_1987_S1640666);


    /// Book histograms and initialise projections before the run
    void init() {
      declare(TriggerUA5(), "Trigger");
      declare(ChargedFinalState(Cuts::abseta < 5.0), "CFS");

      book(_hist_mean_nch, 1, 1, 1);
      book(_hist_nch,      3, 1, 1);
      book(_sumWPassed, "SumW");

      vector<double> tmp; tmp.reserve(41);
      for (int i = 1; i < 41; ++i)  tmp.push_back(2*i-1);
      tmp.push_back(82);
      _axis = YODA::Axis<double>(tmp);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _hist_nch->xEdges();
      // Trigger
      const TriggerUA5& trigger = apply<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;

      _sumWPassed->fill();

      // Count final state particles in several eta regions
      const int Nch = apply<ChargedFinalState>(event, "CFS").size();

      // Fill histograms
      _hist_nch->fill(map2string(Nch));
      _hist_mean_nch->fill(546, Nch);

    }

    string map2string(const size_t nch) const {
      if ( 82 <= nch && nch <=  84)  return   "82.0 - 84.0";
      if ( 86 <= nch && nch <=  88)  return   "86.0 - 88.0";
      if ( 90 <= nch && nch <=  92)  return   "90.0 - 92.0";
      if ( 94 <= nch && nch <=  96)  return   "94.0 - 96.0";
      if ( 98 <= nch && nch <= 100)  return  "98.0 - 100.0";
      if (102 <= nch && nch <= 110)  return "102.0 - 110.0";
      if (112 <= nch && nch <= 120)  return "112.0 - 120.0";
      const size_t idx = _axis.index(nch);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_hist_nch, 1.0 / *_sumWPassed);
      scale(_hist_mean_nch, 1.0 / *_sumWPassed);

      for (auto& b : _hist_nch->bins()) {
        const size_t idx = b.index();
        double sf = idx < 46? 2.0 : 8.0;
        if (idx == 40)  sf = 3.0;
        b.scaleW(1.0/sf);
      }

    }


  private:

    CounterPtr _sumWPassed;

    BinnedHistoPtr<int> _hist_mean_nch;
    BinnedHistoPtr<string> _hist_nch;
    vector<string> _edges;
    YODA::Axis<double> _axis;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(UA5_1987_S1640666, UA5_1987_I244829);

}
