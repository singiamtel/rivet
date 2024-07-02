// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- -> pi+pi-
  class CMD_1985_I221309 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD_1985_I221309);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_sigma1, 1,1,1);
      book(_sigma2, 2,1,1);
      for (const string& en : _sigma1.binning().edges<0>()) {
        double end = std::stod(en)*MeV;
        if(isCompatibleWithSqrtS(end)) {
          ecms = en;
          break;
        }
      }
      icms=-1;
      for(const int& en : _sigma2.binning().edges<0>()) {
        double end = double(en)*MeV;
        if(isCompatibleWithSqrtS(end)) {
          icms = en;
          break;
        }
      }
      if(ecms.empty() && icms<0) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      if(!ecms.empty()) _sigma1->fill(ecms);
      if(icms>0       ) _sigma2->fill(icms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma1,crossSection()/ sumOfWeights() /nanobarn);
      scale(_sigma2,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string>  _sigma1;
    BinnedHistoPtr<int>     _sigma2;
    string ecms;
    int icms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD_1985_I221309);


}
