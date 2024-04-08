// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief f_2, f_0 and K_2 spectra at 29 GeV
  class HRS_1986_I18688 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1986_I18688);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_f2,1,1,1);
      book(_h_f0,1,1,2);
      book(_h_K2,1,1,3);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_f2->xEdges();
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==9010221 ||
                                             Cuts::abspid==225 ||
                                             Cuts::abspid==315)) {
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
      	double xE = 2.*modp/sqrtS();
        if (p.pid()==225) {
          _h_f2->fill(map2string(xE), 1./beta);
        }
        else if (p.pid()==315) {
          _h_K2->fill(map2string(xE), 1./beta);
        }
        else {
          _h_f0->fill(map2string(xE), 1./beta);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale( _h_f0, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale( _h_f2, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale( _h_K2, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      for( auto & hist : {_h_f0,_h_f2,_h_K2} ) {
        for(auto & b: hist->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axis.width(idx));
        }
      }
    }

    ///@}

    string map2string(const double value) const {
      const size_t idx = _axis.index(value);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }

    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_f0, _h_f2, _h_K2;
    YODA::Axis<double> _axis{0.1, 0.2, 0.3, 0.4, 0.7};
    vector<string> _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1986_I18688);

}
