// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief e+ e- > hadrons near Upsilon resonances
  class ARGUS_1994_I375418 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1994_I375418);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(ChargedFinalState(), "CFS");

      // Book histograms
      for(unsigned int ix=0;ix<3;++ix)
        book(_sigma[ix], 1+ix,1,1);
      _eCent = getOption<string>("ECENT", std::to_string(sqrtS()/MeV));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (apply<ChargedFinalState>(event,"CFS").particles().size()<3) vetoEvent;
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // mu+mu- + photons
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        vetoEvent;
      }
      else if (nCount[-11]==1 and nCount[11]==1 && ntotal==2+nCount[22]) {
        vetoEvent;
      }
      else { // everything else
        _sigma[0]->fill(_eCent);
        _sigma[1]->fill(_eCent);
        double H0=0.,H2=0;
        for (unsigned int ix=0; ix<fs.particles().size(); ++ix) {
          double p1 = fs.particles()[ix].p3().mod();
          for (unsigned int iy=0; iy<fs.particles().size(); ++iy) {
            double p2 = fs.particles()[iy].p3().mod();
            double cTheta=fs.particles()[ix].p3().dot(fs.particles()[iy].p3())/p1/p2;
            double pre=p1*p2/s;
            H0 += pre;
            H2 += 0.5*pre*(3.*sqr(cTheta)-1);
          }
        }
        if (H2/H0<0.35) _sigma[2]->fill(_eCent);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/sumOfWeights()/nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _eCent;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1994_I375418);

}
