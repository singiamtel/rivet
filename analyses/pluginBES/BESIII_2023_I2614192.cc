// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- > n nbar
  class BESIII_2023_I2614192 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2614192);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      for(unsigned int ix=0; ix<2; ++ix) {
        book(_c[ix],"/TMP/c_"+toString(ix));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      double cTheta(0.);
      for (const Particle& p :  fs.particles()) {
        if(abs(p.pid())!=PID::NEUTRON) vetoEvent;
        if(p.pid()==PID::NEUTRON) cTheta = p.mom().z()/p.mom().p3().mod();
      }
      _c[0]->fill();
      _c[1]->fill(2.-5.*sqr(cTheta));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // storage of the values to fill histos
      const double mn    = 0.93956542052;
      const double tau   = 4.*sqr(mn/sqrtS());
      const double beta  = sqrt(1.-tau);
      const double alpha = 7.2973525693e-3;
      const double GeV2nb = 0.3893793721e6;
      const double fact = 4./3.*M_PI*sqr(alpha)*beta/2/sqr(sqrtS())*GeV2nb;
      scale(_c, crossSection()/sumOfWeights()/nanobarn/fact);

      double val[3], err[3];
      val[0] = sqrt(2.*_c[1]->val()/(_c[0]->val()-_c[1]->val())/tau);
      err[0] = 0.5/_c[1]->val()/tau*(sqr(_c[0]->err()*_c[1]->val())+sqr(_c[1]->err()*_c[0]->val())) /
               pow(_c[0]->val()-_c[1]->val(),3);
      val[1] = 1e2*sqrt(0.5*(_c[0]->val()-_c[1]->val()));
      err[1] = 1e2*0.125*(sqr(_c[0]->err())+sqr(_c[1]->err()))/(_c[0]->val()-_c[1]->val());
      val[2] = 1e2*sqrt(_c[1]->val()/tau);
      err[2] = 1e2*0.25*sqr(_c[1]->err())/_c[1]->val()/tau;
      // fill histos
      for (unsigned int ix=0; ix<3; ++ix) {
        Estimate1DPtr  mult;
        book(mult, 1, 1, 1+ix);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(val[ix], err[ix]);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2614192);

}
