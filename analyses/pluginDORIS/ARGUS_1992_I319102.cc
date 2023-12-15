// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief charged multiplicity at 4s and nearby continuum
  class ARGUS_1992_I319102 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1992_I319102);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");

      // Book histograms
      if (isCompatibleWithSqrtS(10.47*GeV)) {
        book(_h_N, 2, 1, 1);
        book(_h_tot_N,4,1,1);
      }
      book(_h_N_Upsilon, 3, 1, 1);
      book(_h_N_tot_Upsilon,5,1,1);
      // counters for R
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      book(_w_cont,"/TMP/w_cont");
      book(_w_ups ,"/TMP/w_ups" );
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, unsigned int& nCharged) {
      for (const Particle & p: mother.children()) {
        if (!p.children().empty()) {
          findDecayProducts(p, nCharged);
        }
        else if (PID::isCharged(p.pid())) {
          ++nCharged;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==300553);
      // Continuum
      if (upsilons.empty()) {
        map<long,int> nCount;
        int ntotal(0);
        unsigned int nCharged(0);
        for (const Particle& p : fs.particles()) {
          nCount[p.pid()] += 1;
          ++ntotal;
          if (PID::isCharged(p.pid())) ++nCharged;
        }
        // mu+mu- + photons
        if (nCount[-13]==1 and nCount[13]==1 &&
           ntotal==2+nCount[22])
          _c_muons->fill();
        // everything else
        else {
          _c_hadrons->fill();
          if (_h_N) {
            _h_N->fill(nCharged);
            _h_tot_N->fill(Ecm1, nCharged);
            _w_cont->fill();
          }
        }
      }
      // upsilon 4s
      else {
        for (const Particle& ups : upsilons) {
          unsigned int nCharged(0);
          findDecayProducts(ups,nCharged);
          _h_N_Upsilon->fill(nCharged);
          _h_N_tot_Upsilon->fill(Ecm2, nCharged);
          _w_ups->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate0D R = *_c_hadrons/ *_c_muons;
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Estimate1DPtr hadrons;
      book(hadrons, "sigma_hadrons");
      Estimate1DPtr muons;
      book(muons, "sigma_muons"  );
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 1);
      if (isCompatibleWithSqrtS(9.36*GeV)) {
        mult->bin(1).set(R.val(), R.errPos());
        hadrons->bin(1).set(sig_h, err_h);
        muons  ->bin(1).set(sig_m, err_m);
      }
      if (_h_N) {
      normalize(_h_N, 200.);
      if (_w_cont->val()!=0)
        scale(_h_tot_N,1./ *_w_cont);
      }
      normalize(_h_N_Upsilon, 200.);
      if (_w_ups->val()!=0) {
        scale(_h_N_tot_Upsilon,1./ *_w_ups);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_N, _h_N_Upsilon;
    BinnedHistoPtr<string> _h_tot_N, _h_N_tot_Upsilon;
    CounterPtr _c_hadrons, _c_muons;
    CounterPtr _w_cont,_w_ups;
    const string Ecm1 = "10.47", Ecm2 = "10.575";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1992_I319102);


}
