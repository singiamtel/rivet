// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class BELLE_2005_I677625 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I677625);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if (inRange(sqrtS()/GeV,2.025,4.)) {
        book(_sigmaProton,"TMP/nProton",refData(1,1,1));
        if (inRange(sqrtS()/GeV, 2.075, 2.9) || inRange(sqrtS()/GeV, 3.1, 4.)) {
          double sMin=2.075, sMax=2.1, step=0.1;
          unsigned int ihist=3,iy=1;
          while (sMin<4.) {
            if (inRange(sqrtS()/GeV, sMin, sMax)) break;
            sMin=sMax;
            sMax+=step;
            iy+=1;
            if (iy==4) {
              ihist+=1;
              iy=1;
            }
            if (fuzzyEquals(2.9, sMin)) {
              sMin=3.1;
              sMax=3.5;
              step=0.5;
            }
          }
          book(_h_cTheta[0],ihist,1,iy);
        }
        if (inRange(sqrtS()/GeV, 2.075, 2.5)) {
          book(_h_cTheta[1],2,1,1);
        }
        else if (inRange(sqrtS()/GeV, 2.5, 3.)) {
          book(_h_cTheta[1],2,1,2);
        }
        else if (inRange(sqrtS()/GeV, 3., 4.)) {
          book(_h_cTheta[1],2,1,3);
        }
      }
      else {
        throw Error("Invalid CMS energy for BELLE_2005_I677625");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle & p : part) {
        if (p.pid()==PID::PROTON) {
          foundP=true;
          cTheta = abs(p.mom().z()/p.mom().p3().mod());
        }
        else if (p.pid()==PID::ANTIPROTON) {
          foundM=true;
        }
      }
      if (!foundP || !foundM) vetoEvent;
      if (cTheta<=0.6)    _sigmaProton->fill(sqrtS());
      if (_h_cTheta[0] ) _h_cTheta[0]->fill(cTheta);
      if (_h_cTheta[1] ) _h_cTheta[1]->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      for (unsigned int ix=0; ix<2; ++ix) {
        if (_h_cTheta[ix] ) scale(_h_cTheta[ix], fact);
      }
      scale(_sigmaProton,fact);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_sigmaProton,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cTheta[2];
    Histo1DPtr _sigmaProton;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I677625);

}
