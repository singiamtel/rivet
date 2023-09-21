// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief   F2 structure function in DIS e+ p (ZEUS)

  class ZEUS_1996_I420332 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1996_I420332);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");


      book(_h_f2, {3.2, 4., 5., 7., 9., 11., 13., 16., 20., 32., 40., 50., 65., 85., 110., 140.,
                   185., 240., 310., 410., 530., 710., 900., 1300., 1800., 2500., 3500., 15000.});
      size_t idx = 1;
      for (auto& b : _h_f2->bins()) {
        if (idx == 4 || idx == 7 || idx == 10 || idx == 12 || idx == 14 || idx == 16)  ++idx;
        book(b, idx++, 1, 1);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      //const DISLepton& dl = apply<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;

      // Flux factor
      const double alpha = 7.29927e-3;
      double F = x*sqr(Q2)/(2.0*M_PI*sqr(alpha)*(1.0+sqr(1.-y)));
      _h_f2->fill(Q2,x,F);


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double gev2nb = 0.389e6;
      const double scalefactor=crossSection()/nanobarn/sumOfWeights()/gev2nb ;
      scale(_h_f2, scalefactor);
      divByGroupWidth(_h_f2);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_f2;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ZEUS_1996_I420332);

}
