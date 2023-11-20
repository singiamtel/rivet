// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief rho0, K*=, K*0 spectra at 29 GeV
  class HRS_1989_I276948 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1989_I276948);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h["rho"], 1, 1, 1);
      book(_h["Kstar0"], 2, 1, 1);
      book(_h["Kstarp"], 2, 1, 2);
      book(_h_sig_rho    ,3,1,1);
      book(_h_sig_Kstar0 ,4,1,1);
      book(_h_mult_rho   ,3,1,3);
      book(_h_mult_Kstar0,4,1,3);
      _axes["rho"] = YODA::Axis<double>({0.0725, 0.0925, 0.114, 0.136, 0.159, 0.183, 0.207, 0.231,
                                         0.2555, 0.28, 0.3045, 0.3295, 0.354, 0.3785, 0.4035, 0.4285,
                                         0.4595, 0.5155, 0.6025, 0.702, 0.802, 0.902});
      _axes["Kstar0"] = YODA::Axis<double>({0.077, 0.119, 0.163, 0.2095, 0.2575, 0.3065,
                                            0.3555, 0.405, 0.455, 0.5, 0.6, 0.7, 1.0});
      _axes["Kstarp"] = _axes["Kstar0"];
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty()) {
        for (const auto& item : _h) {
          _edges[item.first] = item.second->xEdges();
        }
      }
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==113 || Cuts::abspid==323 || Cuts::abspid==313)) {
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
      	double xE = 2.*modp/sqrtS();
        if (p.pid()==113) {
          discfill("rho", xE, 1./beta);
          _h_sig_rho->fill(29);
          _h_mult_rho->fill(29);
        }
        else if (p.pid()==313) {
          discfill("Kstar0", xE, 1./beta);
          _h_sig_Kstar0->fill(29);
          _h_mult_Kstar0->fill(29);
        }
        else {
          discfill("Kstarp", xE, 1./beta);
        }
      }
    }

    void discfill(const string& name, const double val, const double weight) {
      string edge = "OTHER";
      const size_t idx = _axes[name].index(val);
      if (idx && idx <= _edges[name].size())  edge = _edges[name][idx-1];
      _h[name]->fill(edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale({_h_sig_rho, _h_sig_Kstar0}, crossSection()/picobarn/sumOfWeights());
      scale({_h_mult_rho, _h_mult_Kstar0}, 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<int> _h_sig_rho,_h_sig_Kstar0;
    BinnedHistoPtr<int> _h_mult_rho,_h_mult_Kstar0;
    map<string,BinnedHistoPtr<string>> _h;
    map<string,YODA::Axis<double>> _axes;
    map<string,vector<string>> _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1989_I276948);

}
