// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
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
      vector<string> energies({"400.0", "440.0", "642.6", "662.6", "677.6", "687.6", "697.6", "707.6", "717.6",
          "727.6", "737.6", "747.6", "757.6", "767.6", "773.6", "775.6", "777.6", "779.6", "781.6", "783.6",
          "785.6", "787.6", "789.6", "791.6", "797.6", "807.6", "817.6", "827.6", "837.6", "847.6", "857.6",
          "867.6", "877.6", "887.6", "897.6", "907.6", "917.6", "927.6", "937.6", "947.6", "957.6", "967.6",
          "977.6", "987.6", "997.6", "1006.6", "1067.0", "1077.0", "1087.0", "1097.0", "1107.0", "1117.0",
          "1127.0", "1137.0", "1152.0", "1167.0", "1177.0", "1187.0", "1197.0", "1207.0", "1217.0", "1227.0",
          "1237.0", "1247.0", "1257.0", "1267.0", "1277.0", "1287.0", "1297.0", "1307.0", "1317.0", "1327.0",
          "1337.0", "1347.0", "1357.0", "1367.0", "1377.0", "1387.0", "1397.0"});
      for(const string& en : energies) {
        double end = std::stod(en)*MeV;
        if(isCompatibleWithSqrtS(end)) {
          ecms = en;
          break;
        }
      }
      vector<int> enint({360, 380, 410, 430, 438, 470, 540, 580, 620, 640, 660,
          700, 740, 760, 770, 774, 778, 780, 782, 786, 790, 794, 800, 820});
      icms=-1;
      for(const int& en : enint) {
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
