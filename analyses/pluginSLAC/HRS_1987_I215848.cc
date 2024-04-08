// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Hadron Spectra in $e^+e^-$ collisions at 29 GeV
  class HRS_1987_I215848 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1987_I215848);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // hists
      book(_d["pi"],  2, 1, 1);
      book(_d["Kp"],  3, 1, 1);
      book(_d["p"],   4, 1, 1);
      book(_h["K0"],  5, 1, 1);
      book(_h["lam"], 6, 1, 1);
      _axes["pi"]=YODA::Axis<double>({0.035800994999999995,0.049226100999999994,
          0.06281088,0.076470283,0.090170398,0.103895139,0.11763585,0.13138755});
      _axes["Kp"]=YODA::Axis<double>({0.048458616,0.059073984,0.070793603,0.083151849,
          0.095902214,0.108906933,0.122084929,0.135385471});
      _axes["p"]=YODA::Axis<double>({0.073322849,0.08073250800000001,0.08966457400000001,
          0.099710693,0.110567797,0.122019328,0.133913114,0.14614108599999998});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty()) {
        for (const auto& item : _d) {
          _edges[item.first] = item.second->xEdges();
        }
      }
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle & p : ufs.particles()) {
        const double xE = 2.*p.E()/sqrtS();
        const double beta = p.p3().mod() / p.E();
        if (p.pid()==130 || p.pid()==310) {
          _h["K0"]->fill(xE,1./beta);
        }
        else if (p.abspid()==321) {
          discfill("Kp",xE, 1./beta);
        }
        else if (p.abspid()==211) {
          discfill("pi",xE, 1./beta);
        }
        else if (p.abspid()==2212) {
          discfill("p",xE, 1./beta);
        }
        else if (p.abspid()==3122) {
          _h["lam"]->fill(xE, 1./beta);
        }
      }
    }

    void discfill(const string& name, const double value, const double weight) {
      string edge = "OTHER";
      const size_t idx = _axes[name].index(value);
      if (idx && idx <= _edges[name].size())  edge = _edges[name][idx-1];
      _d[name]->fill(edge,weight);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()*sqr(sqrtS())/microbarn/sumOfWeights();
      scale(_h, sf);
      scale(_d, sf);
      for( auto & hist : _d) {
        for(auto & b: hist.second->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axes[hist.first].width(idx));
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, BinnedHistoPtr<string> > _d;
    map<string, YODA::Axis<double> > _axes;
    map<string, vector<string> > _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1987_I215848);

}
