// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi production at 13 TeV (central rapidity)
  class ALICE_2021_I1898832 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1898832);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(Cuts::pid==443), "UFS");
      book(_h_pT,1,1,1);
      book(_h_y1,2,1,1);
      book(_h_y2,3,1,1);
      _axis = YODA::Axis<double>({ 0., 0.5, 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 30., 40.});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_pT->xEdges();
      // loop over J/Psi
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        const double absrap = p.absrap();
        if (absrap>0.9) continue;

        _h_pT->fill(disc(p.perp()));

        _h_y1->fill("0.0"s);

        string yedge("0.7");
        if (absrap < 0.2) yedge = "0.1";
        else if (absrap < 0.5) yedge = "0.3";

        _h_y2->fill(yedge);
      }
    }

    string disc(const double value) const {
      string edge("OTHER");
      const size_t idx = _axis.index(value);
      if (0 < idx && idx <= _axis.numBins())  return _edges[idx-1];
      return edge;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/microbarn/sumOfWeights();
      // 1.8 for the rapidity -0.9 to 0.9 in the double differential
      scale(_h_pT,fact/1.8);
      for(auto & b : _h_pT->bins()) {
        b.scaleW(1./_axis.width(b.index()));
      }
      // divide by range -0.9 to 0.9 as quoted as dsigma/dy
      scale(_h_y1,fact/1.8);
      // 0.5 as fold + and - rapidity
      scale(_h_y2,0.5*fact);
      for(auto & b : _h_y2->bins()) {
        if     (b.index()==1) b.scaleW(1./0.2);
        else if(b.index()==2) b.scaleW(1./0.3);
        else                  b.scaleW(1./0.4);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_pT,_h_y1,_h_y2;
    vector<string> _edges;
    YODA::Axis<double> _axis;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1898832);

}
