// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi0 pi0
  class BELLE_2009_I815978 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I815978);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if (inRange(sqrtS()/GeV,0.6,3.3) || inRange(sqrtS()/GeV,3.6,4.1)) {
        if (sqrtS()>0.72) book(_sigmapipi[0],"TMP/npipi_1",refData(31, 1, 1));
        book(_sigmapipi[1],"TMP/npipi_2",refData(31, 1, 2));
        double sMin=0.6, step=0.02;
        unsigned int ihist=1,iy=1;
        while (sMin<4.1) {
          if (inRange(sqrtS()/GeV, sMin, sMin+step)) {
            break;
          }
          sMin+=step;
          iy+=1;
          if (iy==4) {
            ihist+=1;
            iy=1;
          }
          if(fuzzyEquals(1.8, sMin)) step=0.04;
          else if(fuzzyEquals(2.4, sMin)) step=0.1;
          else if(fuzzyEquals(3.2, sMin)) sMin=3.6;
        }
        if (!inRange(sqrtS()/GeV,3.2,3.6)) book(_h_cTheta,ihist,1,iy);
      }
      else {
       	throw Error("Invalid CMS energy for BELLE_2009_I815978");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      for (const Particle & p : part) {
        if (p.pid()!=PID::PI0) vetoEvent;
      }
      double cTheta = abs(part[0].mom().z()/part[0].mom().p3().mod());
      if (cTheta<=0.6&&_sigmapipi[0]) _sigmapipi[0]->fill(sqrtS());
      if (cTheta<=0.8&&_sigmapipi[1]) _sigmapipi[1]->fill(sqrtS());
      if (_h_cTheta ) _h_cTheta ->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) scale(_h_cTheta ,fact);
      for (unsigned int ix=0; ix<2; ++ix) {
        if (!_sigmapipi[ix]) continue;
        scale(_sigmapipi[ix],fact);
        Estimate1DPtr tmp;
        book(tmp, 31, 1, 1+ix);
        barchart(_sigmapipi[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmapipi[2],_h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I815978);

}
