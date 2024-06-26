// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged particle spectra between 2.2 and 4.8 GeV
  class BESII_2004_I622224 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2004_I622224);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState fs;
      declare(fs, "FS");
      unsigned int iloc(0);
      if(isCompatibleWithSqrtS(2.2 , 1E-3)) iloc = 1;
      else if(isCompatibleWithSqrtS(2.6 , 1E-3)) iloc = 2;
      else if(isCompatibleWithSqrtS(3.0 , 1E-3)) iloc = 3;
      else if(isCompatibleWithSqrtS(3.2 , 1E-3)) iloc = 4;
      else if(isCompatibleWithSqrtS(4.6 , 1E-3)) iloc = 5;
      else if(isCompatibleWithSqrtS(4.8 , 1E-3)) iloc = 6;
      else MSG_ERROR("Beam energy not supported!");
      assert(iloc != 0);
      book(_h_ln, iloc, 1, 1);
      book(_h_weight, "TMP/Weight");
    }

    string map2string(const double xi) const {
      const size_t idx = axis.index(xi) - 1;
      if (idx < edges.size())  return edges[idx];
      return "OTHER";
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (edges.empty())  edges = _h_ln->xEdges();
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      if (fs.particles().size()==2 &&
          abs(fs.particles()[0].pid())==13 &&
          abs(fs.particles()[1].pid())==13) vetoEvent;
      for (const Particle& p : fs.particles()) {
        const Vector3 mom3 = p.p3();
        double pp = mom3.mod();
        double xi = -log(2.*pp/sqrtS());
        _h_ln->fill(map2string(xi));
      }
      _h_weight->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_ln, 1./dbl(*_h_weight));
      for(auto & b: _h_ln->bins()) {
        const size_t idx = b.index();
        b.scaleW(1./axis.width(idx));
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_ln;
    CounterPtr _h_weight;
    vector<string> edges;
    YODA::Axis<double> axis{0.05, 0.45, 0.6, 0.7, 0.8, 0.9, 1.0,
                            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                            2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2004_I622224);

}
