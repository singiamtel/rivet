// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > pi+ pi-
  class CMD3_2023_I2634277 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2023_I2634277);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_sigma[ix], 1, 1+ix, 1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*MeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty() && _ecms[2].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if (abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      for(unsigned int ix=0;ix<3;++ix) {
        if(!_ecms[ix].empty()) _sigma[ix]->fill(_ecms[ix]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double gev2nb =0.389e6;
      const double beta = sqrt(1.-4.*sqr(0.13957039/sqrtS()));
      const double alpha = 1./137.035999084;
      const double sigma0 = gev2nb*M_PI*sqr(alpha)*beta*sqr(beta)/3/sqr(sqrtS());
      for(unsigned int ix=0;ix<3;++ix)
        scale(_sigma[ix], crossSection()/ sumOfWeights() /nanobarn /sigma0);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2023_I2634277);

}
