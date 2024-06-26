// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> K0S K0S
  class BELLE_2013_I1245023 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1245023);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if (inRange(sqrtS()/GeV,1.05,4.0)) {
        if (sqrtS()>1.1) book(_sigmaKK[0],"TMP/nKK_1",refData( 1, 1, 1));
        book(_sigmaKK[1],"TMP/nKK_2",refData( 1, 1, 2));
        double sMin=1.1, step=0.01;
        unsigned int ihist=2,iy=1;
        while (sMin<3.3) {
          if (inRange(sqrtS()/GeV, sMin, sMin+step)) {
            break;
          }
          sMin+=step;
          iy+=1;
          if (iy==4) {
            ihist+=1;
            iy=1;
          }
          if (fuzzyEquals(1.9, sMin)) step=0.02;
          else if (fuzzyEquals(2.4, sMin)) step=0.04;
          else if (fuzzyEquals(2.6, sMin)) step=0.1;
        }
        if (ihist<=40) book(_h_cTheta,ihist,1,iy);
      }
      else {
       	throw Error("Invalid CMS energy for BELLE_2013_I1245023");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      for (const Particle & p : part) {
        if (p.pid()!=PID::K0S) vetoEvent;
      }
      const double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if (cTheta<=0.6)               _sigmaKK[1]->fill(sqrtS());
      if (cTheta<=0.8 &&_sigmaKK[0]) _sigmaKK[0]->fill(sqrtS());
      if (_h_cTheta )                _h_cTheta ->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) scale(_h_cTheta ,fact);
      for(unsigned int ix=0; ix<2; ++ix) {
        if (!_sigmaKK[ix]) continue;
        scale(_sigmaKK[ix],fact);
        Estimate1DPtr tmp;
        book(tmp,1,1,1+ix);
        barchart(_sigmaKK[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaKK[2],_h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1245023);

}
