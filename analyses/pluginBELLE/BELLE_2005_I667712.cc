// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// gamma gamma -> pi+pi-/K+ K- for energies between 2.4 and 4.1 GeV
  class BELLE_2005_I667712 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I667712);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if (sqrtS()<2.4*GeV || sqrtS()>4.1*GeV)
        throw Error("Invalid CMS energy for BELLE_2005_I667712");
      for (unsigned int ix=0;ix<7;++ix) {
        std::ostringstream title;
        title << "/TMP/nPi_" << ix;
        book(_cPi[ix], title.str());
      }
      for (unsigned int ix=0;ix<7;++ix) {
        std::ostringstream title;
        title << "/TMP/nK_" << ix;
        book(_cK[ix], title.str());
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size() != 2) vetoEvent;
      if (part[0].pid() != -part[1].pid()) vetoEvent;
      double cTheta(0.);
      bool foundPi(false),foundK(false);
      for (const Particle& p : part) {
        if (p.pid() == PID::PIPLUS) {
          foundPi=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid() == PID::KPLUS) {
          foundK=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
      }
      if (!foundPi && !foundK) vetoEvent;
      int ibin = cTheta/0.1;
      if (ibin>5) vetoEvent;
      if (foundPi) {
        _cPi[   0  ]->fill();
        _cPi[ibin+1]->fill();
      }
      else if (foundK) {
        _cK[   0  ]->fill();
        _cK[ibin+1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for (unsigned int ip=0;ip<2;++ip) {
        CounterPtr denom = ip==0 ? _cPi[0] : _cK[0];
        if (denom->numEntries()==0) continue;
        for (unsigned int ih=0;ih<7;++ih) {
          CounterPtr numer = ip==0 ? _cPi[ih] : _cK[ih];
          double sigma = numer->val()*fact;
          double error = numer->err()*fact;
          // bin width for 2d dist
          if (ih!=0) {
            sigma /=0.1;
            error /=0.1;
          }
          unsigned int ix=5+ip, iy=ih;
          if (ih==0) {
            ix=1;
            iy = ip==0 ? 2 : 1;
          }
          // ratio
          Estimate0D temp;
          Estimate1DPtr cross, ratio;
          book(cross, ix, 1, iy);
          if (ih!=0) {
            book(ratio,ix-2,1,iy);
            temp = *numer / *denom;
          }
          for (auto& b : cross->bins()) {
            if (inRange(sqrtS(), b.xMin(), b.xMax())) {
              b.set(sigma, error);
              if (ih!=0) {
                ratio->bin(b.index()).set(temp.val()/0.1, temp.errPos()/0.1);
              }
            }
          }
        }
      }

    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _cPi[7],_cK[7];
    ///@}

  };



  RIVET_DECLARE_PLUGIN(BELLE_2005_I667712);

}
