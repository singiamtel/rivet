// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief  gamma, pi0 and eta spectra at 14, 22.5 and 34.4 GeV
  class JADE_1985_I213948 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1985_I213948);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // find the beam energy
      int ioff=-1;
      if (isCompatibleWithSqrtS(34.5*GeV))  ioff=0;
      else if(isCompatibleWithSqrtS(22.5*GeV)) 	ioff=1;
      else if(isCompatibleWithSqrtS(14.0*GeV)) 	ioff=2;
      else {
      	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");
      }
      // book histos
      book(_h["gamma"], ioff+1, 1, 1);
      book(_h["pi0"],   ioff+4, 1, 1);
      if (ioff==0) {
        book(_h["eta"], 7, 1, 1);
        _axes.emplace_back<vector<double>>({0.00, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                                            0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20,
                                            0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40,
                                            0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90}); // d01
        _axes.emplace_back<vector<double>>({0.011, 0.023, 0.035, 0.047, 0.057, 0.069, 0.079, 0.083,
                                            0.093, 0.107, 0.113, 0.119, 0.143, 0.167, 0.191, 0.239}); // d04
        _axes.emplace_back<vector<double>>({0.040, 0.056, 0.088, 0.274}); // d07
      }
      else if (ioff==2) {
        _axes.emplace_back<vector<double>>({0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                                            0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.40, 0.60, 1.00}); // d02
        _axes.emplace_back<vector<double>>({0.018, 0.036, 0.054, 0.072, 0.108, 0.144, 0.216, 0.360}); // d05
      }
      else {
        _axes.emplace_back<vector<double>>({0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                                            0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.40, 0.60, 1.00}); //d03
        _axes.emplace_back<vector<double>>({0.030, 0.058, 0.086, 0.114, 0.170, 0.226, 0.338, 0.562}); // d06
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty()) {
        _edges.resize(_axes.size());
        _edges[0] = _h["gamma"]->xEdges();
        _edges[1] = _h["pi0"]->xEdges();
        if (_h.count("eta"))  _edges[3] = _h["eta"]->xEdges();
      }
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = 0.5*(beams.first.p3().mod()+beams.second.p3().mod());
      // gamma
      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles(Cuts::pid==22)) {
      	const double xE = p.E()/meanBeamMom;
        _h["gamma"]->fill(map2string(xE, 0));
      }
      // pi0, eta
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==111 or Cuts::pid==221)) {
      	double xE = p.E()/meanBeamMom;
      	if (p.pid()==111) {
          _h["pi0"]->fill(map2string(xE, 1));
        }
      	else if(_h.count("eta")) {
          _h["eta"]->fill(map2string(xE, 2));
        }
      }
    }

    string map2string(const double val, const size_t k) const {
      const size_t idx = _axes[k].index(val);
      if (idx && idx <= _axes[k].numBins())  return _edges[k][idx-1];
      return "OTHER";
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()*sqr(sqrtS())/microbarn/sumOfWeights();
      scale(_h, sf);
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, BinnedHistoPtr<string>> _h;
    vector<YODA::Axis<double>> _axes;
    vector<vector<string>> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(JADE_1985_I213948);

}
