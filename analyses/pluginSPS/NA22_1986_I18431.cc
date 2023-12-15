// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// NA22 min bias multiplicity distributions.
  class NA22_1986_I18431 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NA22_1986_I18431);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(ChargedFinalState(), "CFS");

      // Figure out beam type
      const ParticlePair& beam = beams();
      int btype = 0;
      if (beam.first.pid() == PID::PIPLUS && beam.second.pid() == PID::PROTON) {
        btype = 1;
        sf *= 20.94; // simga_inel(pi+ p) from DOI:10.1007/BF01550769
      }
      else if (beam.first.pid() == PID::KPLUS && beam.second.pid() == PID::PROTON) {
        btype = 2;
        sf *= 17.72; // simga_inel(K+ p) from DOI:10.1007/BF01550769
      }
      else if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON) {
        btype = 3;
        sf *= 32.40; // simga_inel(p p) from DOI:10.1007/BF01550769
      }
      else {
        MSG_ERROR("Beam error: Not compatible!");
        return;
      }

      // Book histo for appropriate beam type
      if (btype == 3)  book(_h_mult2, btype, 1, 1);
      else             book(_h_mult1, btype, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      size_t nfs = apply<ChargedFinalState>(event, "CFS").size();
      if (nfs >= 30)  nfs = 30;
      if (_h_mult1) {
        string edge = to_string(nfs)+".0";
        if (nfs == 30)  edge = ">= 30.0";
        _h_mult1->fill(edge);
      }
      if (_h_mult2)  _h_mult2->fill(nfs);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_h_mult1) normalize(_h_mult1, sf); // normalize to 2 to account for bin width
      if (_h_mult2) normalize(_h_mult2, sf); // normalize to 2 to account for bin width
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_mult1;
    BinnedHistoPtr<int> _h_mult2;
    double sf = 2.0;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(NA22_1986_I18431);


}
