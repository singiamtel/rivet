// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  /// UA5 min bias charged multiplicities in central \f$ \eta \f$ ranges
  class UA5_1989_I267179 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(UA5_1989_I267179);


    /// @name Analysis methods
    /// @{

    /// Book histograms and projections
    void init() {
      declare(TriggerUA5(), "Trigger");
      declare(ChargedFinalState(Cuts::abseta < 0.5), "CFS05");
      declare(ChargedFinalState(Cuts::abseta < 1.5), "CFS15");
      declare(ChargedFinalState(Cuts::abseta < 3.0), "CFS30");
      declare(ChargedFinalState(Cuts::abseta < 5.0), "CFS50");

      // NB. _hist_nch and _hist_ncheta50 use the same data but different binning
      vector<double> tmp;
      if (isCompatibleWithSqrtS(200*GeV)) {
        book(_hist_nch,        1, 1, 1);
        book(_hist_nch_eta05,  3, 1, 1);
        book(_hist_nch_eta15,  4, 1, 1);
        book(_hist_nch_eta30,  5, 1, 1);
        book(_hist_nch_eta50,  6, 1, 1);
        book(_hist_mean_nch , 11, 1, 1);
        Ecm = 200;

        // _hist_nch axis
        tmp.clear(); tmp.reserve(32);
        for (size_t i = 1; i < 30; ++i)  tmp.push_back(2*i-1);
        tmp.push_back(60); tmp.push_back(62); tmp.push_back(74);
        _axes[0] = YODA::Axis<double>(tmp);

        // _hist_nch_eta05 axis
        _axes[1] = YODA::Axis<double>({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,
                                       7.5, 8.5, 9.5, 10.5, 11.5, 13.0, 19.0});

        // _hist_nch_eta15 axis
        tmp.clear(); tmp.reserve(33);
        for (size_t i = 1; i < 31; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(31); tmp.push_back(32); tmp.push_back(45);
        _axes[2] = YODA::Axis<double>(tmp);

        // _hist_nch_eta30 axis
        tmp.clear(); tmp.reserve(54);
        for (size_t i = 0; i < 49; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(49); tmp.push_back(50);
        tmp.push_back(51); tmp.push_back(53);
        tmp.push_back(81);
        _axes[3] = YODA::Axis<double>(tmp);

        // _hist_nch_eta50 axis
        tmp.clear(); tmp.reserve(59);
        tmp.push_back(0);
        for (size_t i = 3; i < 56; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(56); tmp.push_back(57);
        tmp.push_back(58); tmp.push_back(60);
        tmp.push_back(74);
        _axes[4] = YODA::Axis<double>(tmp);
      }
      else if (isCompatibleWithSqrtS(900*GeV)) {
        book(_hist_nch,        2, 1, 1);
        book(_hist_nch_eta05,  7, 1, 1);
        book(_hist_nch_eta15,  8, 1, 1);
        book(_hist_nch_eta30,  9, 1, 1);
        book(_hist_nch_eta50, 10, 1, 1);
        book(_hist_mean_nch,  12, 1, 1);
        Ecm = 900;

        // _hist_nch axis
        tmp.clear(); tmp.reserve(55);
        for (size_t i = 4; i < 53; ++i)  tmp.push_back(2*i-1);
        tmp.push_back(106); tmp.push_back(108); tmp.push_back(110);
        tmp.push_back(112); tmp.push_back(116);
        _axes[0] = YODA::Axis<double>(tmp);

        // _hist_nch_eta05 axis
        tmp.clear(); tmp.reserve(23);
        for (size_t i = 0; i < 21; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(21); tmp.push_back(23); tmp.push_back(30);
        _axes[1] = YODA::Axis<double>(tmp);

        // _hist_nch_eta15 axis
        tmp.clear(); tmp.reserve(53);
        for (size_t i = 0; i < 47; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(47); tmp.push_back(48); tmp.push_back(49);
        tmp.push_back(51); tmp.push_back(55); tmp.push_back(71);
        _axes[2] = YODA::Axis<double>(tmp);

        // _hist_nch_eta30 axis
        tmp.clear(); tmp.reserve(84);
        for (size_t i = 0; i < 76; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(76); tmp.push_back(77); tmp.push_back(78);
        tmp.push_back(79); tmp.push_back(81); tmp.push_back(84);
        tmp.push_back(90); tmp.push_back(113);
        _axes[3] = YODA::Axis<double>(tmp);

        // _hist_nch_eta50 axis
        tmp.clear(); tmp.reserve(102);
        tmp.push_back(0); tmp.push_back(2);
        for (size_t i = 4; i < 95; ++i)  tmp.push_back(i-0.5);
        tmp.push_back(95);  tmp.push_back(96);  tmp.push_back(97);
        tmp.push_back(98);  tmp.push_back(99);  tmp.push_back(102);
        tmp.push_back(106); tmp.push_back(112); tmp.push_back(130);
        _axes[4] = YODA::Axis<double>(tmp);
      }
      book(_sumWPassed, "SumW");
      /// @todo Moments of distributions
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = apply<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;

      _sumWPassed->fill();

      // Count final state particles in several eta regions
      const int numP05 = apply<ChargedFinalState>(event, "CFS05").size();
      const int numP15 = apply<ChargedFinalState>(event, "CFS15").size();
      const int numP30 = apply<ChargedFinalState>(event, "CFS30").size();
      const int numP50 = apply<ChargedFinalState>(event, "CFS50").size();

      // Fill histograms
      _hist_nch->fill(map2string(numP50, 0));
      _hist_nch_eta05->fill(map2string(numP05, 1));
      _hist_nch_eta15->fill(map2string(numP15, 2));
      _hist_nch_eta30->fill(map2string(numP30, 3));
      _hist_nch_eta50->fill(map2string(numP50, 4));
      _hist_mean_nch->fill(Ecm, numP50);
    }

    string map2string(const int n, const int b) const {
      if (Ecm == 200) {
        if (b == 0) {
          if (64 <= n && n <= 76)  return "64.0 - 76.0";
          if (60 <= n && n <= 62)  return "60.0 - 62.0";
          if (58 == n)             return "58";
        }
        if (b == 1) {
          if (13 <= n && n <= 19)  return "13.0 - 19.0";
        }
        if (b == 2) {
          if (33 <= n && n <= 45)  return "33.0 - 45.0";
          if (31 <= n && n <= 32)  return "31.0 - 32.0";
          if (30 == n)             return "30.0";
        }
        if (b == 3) {
          if (56 <= n && n <= 74)  return "56.0 - 74.0";
          if (53 <= n && n <= 55)  return "53.0 - 55.0";
          if (51 <= n && n <= 52)  return "51.0 - 52.0";
          if (49 <= n && n <= 50)  return "49.0 - 50.0";
          if (48 == n)             return "48.0";
        }
        if (b == 4) {
          if (63 <= n && n <= 77)  return "63.0 - 77.0";
          if (60 <= n && n <= 62)  return "60.0 - 62.0";
          if (58 <= n && n <= 59)  return "58.0 - 59.0";
          if (56 <= n && n <= 57)  return "56.0 - 57.0";
          if (55 == n)             return "55.0";
          if ( 0 <= n && n <=  2)  return "0.0 - 2.0";
        }
      }
      if (Ecm == 600) {
        if (b == 0) {
          if (122 <= n && n <= 136)  return "122.0 - 136.0";
          if (114 <= n && n <= 120)  return "114.0 - 120.0";
          if (106 <= n && n <= 108)  return "106.0 - 108.0";
          if (104 == n)              return "104.0";
          if (  2 <= n && n <=   6)  return "2.0 - 6.0";
        }
        if (b == 1) {
          if (24 <= n && n <= 31)  return "24.0 - 31.0";
          if (21 <= n && n <= 23)  return "21.0 - 23.0";
          if (20 == n)             return "20.0";
        }
        if (b == 2) {
          if (59 <= n && n <= 75)  return "59.0 - 75.0";
          if (54 <= n && n <= 58)  return "54.0 - 58.0";
          if (51 <= n && n <= 53)  return "51.0 - 53.0";
          if (49 <= n && n <= 50)  return "49.0 - 50.0";
          if (47 <= n && n <= 48)  return "47.0 - 48.0";
          if (46 == n)             return "46.0";
        }
        if (b == 3) {
          if (96 <= n && n <= 119)  return "69.0 - 119.0";
          if (89 <= n && n <=  95)  return "89.0 - 95.0";
          if (85 <= n && n <=  88)  return "85.0 - 88.0";
          if (82 <= n && n <=  84)  return "82.0 - 84.0";
          if (80 <= n && n <=  81)  return "80.0 - 81.0";
          if (78 <= n && n <=  79)  return "78.0 - 79.0";
          if (76 <= n && n <=  77)  return "76.0 - 77.0";
          if (75 == n)              return "75.0";
        }
        if (b == 4) {
          if (119 <= n && n <=  137)  return "119.0 - 137.0";
          if (112 <= n && n <=  118)  return "112.0 - 118.0";
          if (107 <= n && n <=  111)  return "107.0 - 111.0";
          if (103 <= n && n <=  106)  return "103.0 - 106.0";
          if (101 <= n && n <=  102)  return "101.0 - 102.0";
          if ( 99 <= n && n <=  100)  return "99.0 - 100.0";
          if ( 97 <= n && n <=   98)  return "97.0 - 98.0";
          if ( 95 <= n && n <=   96)  return "95.0 - 96.0";
          if ( 94 == n)               return "94.0";
          if (  0 <= n && n <=    2)  return "0.0 - 2.0";
        }
      }
      const size_t idx = _axes[b].index((double)n);
      double val = _axes[b].mid(idx);
      val = val + 0.5 - (val<0);
      return std::to_string((int)val) + ".0";
    }

    void finalize() {
      scale(_hist_nch, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta05, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta15, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta30, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta50, 1.0 / *_sumWPassed);
      scale(_hist_mean_nch,  1.0 / *_sumWPassed);
      for (auto& b: _hist_nch->bins()) {
        const size_t idx = b.index();
        double sf = _axes[0].width(idx);
        b.scaleW(1.0/sf);
      }
      for (auto& b: _hist_nch_eta05->bins()) {
        const size_t idx = b.index();
        double sf = _axes[1].width(idx);
        b.scaleW(1.0/sf);
      }
      for (auto& b: _hist_nch_eta15->bins()) {
        const size_t idx = b.index();
        double sf = _axes[2].width(idx);
        b.scaleW(1.0/sf);
      }
      for (auto& b: _hist_nch_eta30->bins()) {
        const size_t idx = b.index();
        double sf = _axes[3].width(idx);
        b.scaleW(1.0/sf);
      }
      for (auto& b: _hist_nch_eta50->bins()) {
        const size_t idx = b.index();
        double sf = _axes[4].width(idx);
        b.scaleW(1.0/sf);
      }
    }

    /// @}


  private:

    /// Weight counter
    CounterPtr _sumWPassed;

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _hist_nch;
    BinnedHistoPtr<string> _hist_nch_eta05;
    BinnedHistoPtr<string> _hist_nch_eta15;
    BinnedHistoPtr<string> _hist_nch_eta30;
    BinnedHistoPtr<string> _hist_nch_eta50;
    BinnedHistoPtr<int> _hist_mean_nch;
    YODA::Axis<double> _axes[5];
    int Ecm;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(UA5_1989_I267179, UA5_1989_S1926373);

}
