// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class LENA_1982_I179431 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LENA_1982_I179431);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for (size_t i=0; i<2; ++i) {
        const auto& ref = refData(i+1, 1, 1);
        book(_h_hadrons[i], "d0"+to_string(i+1)+"_sigma_hadrons", ref);
        book(_h_muons[i], "d0"+to_string(i+1)+"_sigma_muons", ref);
        book(_ratios[i], i+1, 1, 1);
      }
      book(_c_hadrons, "d03_sigma_hadrons");
      book(_c_muons, "d03_sigma_muons");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // mu+mu- + photons
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        _h_muons[0]->fill(sqrtS()/GeV);
        _h_muons[1]->fill(sqrtS()/GeV);
        _c_muons->fill();
      }
      // everything else
      else {
        _h_hadrons[0]->fill(sqrtS()/GeV);
        _h_hadrons[1]->fill(sqrtS()/GeV);
        _c_hadrons->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection() / sumOfWeights() / picobarn;
      scale(_h_hadrons, fact);
      scale(_h_muons, fact);
      scale({_c_hadrons, _c_muons}, fact);
      divide(_h_hadrons[0], _h_muons[0], _ratios[0]);
      divide(_h_hadrons[1], _h_muons[1], _ratios[1]);
      const Estimate0D R =  *_c_hadrons / *_c_muons;
      BinnedEstimatePtr<string> mult;
      book(mult, 3, 1, 1);
      for (auto& b : mult->bins()) {
        const double Ecm = std::stod(b.xEdge());
        if (isCompatibleWithSqrtS(Ecm/GeV)) {
          b.set(R.val(), R.errPos());
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons, _c_muons;
    Histo1DPtr _h_hadrons[2], _h_muons[2];
    Estimate1DPtr _ratios[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LENA_1982_I179431);


}
