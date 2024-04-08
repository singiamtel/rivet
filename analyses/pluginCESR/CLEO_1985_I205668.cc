// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Spectra in Upsilon(1S) decay and nearby continuum
  class CLEO_1985_I205668 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1985_I205668);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      book(_weightSum_cont,"TMP/weightSumcont");
      book(_weightSum_Ups1,"TMP/weightSumUps1");
      // multiplcities
      for (size_t ix=0; ix<2; ++ix) {
        for (size_t iy=0; iy<12; ++iy) {
          book(_mult[ix][iy],"/TMP/MULT_" +toString(ix) + "_" +toString(iy));
        }
      }
      // cont spectra
      book(_cont["pip"]   , 1,1,1);
      book(_cont["Kp"]    , 2,1,1);
      book(_cont["p"]     , 3,1,1);
      book(_cont["pi0"]   , 4,1,1);
      book(_cont["K0"]    , 5,1,1);
      book(_cont["lam"]   , 6,1,1);
      book(_cont["xi"]    , 7,1,1);
      book(_cont["rho"]   , 8,1,1);
      book(_cont["Kstarp"], 9,1,1);
      book(_cont["Kstar0"],10,1,1);
      book(_cont["phi"]   ,11,1,1);
      // ups spectra
      book(_ups1["pip"]   , 1,1,2);
      book(_ups1["Kp"]    , 2,1,2);
      book(_ups1["p"]     , 3,1,2);
      book(_ups1["pi0"]   , 4,1,2);
      book(_ups1["K0"]    , 5,1,2);
      book(_ups1["lam"]   , 6,1,2);
      book(_ups1["xi"]    , 7,1,2);
      book(_ups1["rho"]   , 8,1,2);
      book(_ups1["Kstarp"], 9,1,2);
      book(_ups1["Kstar0"],10,1,2);
      book(_ups1["phi"]   ,11,1,2);

      _axes[0]["pip"] = YODA::Axis<double>({0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17,
                                            0.19, 0.58, 0.68, 0.78, 0.98});
      _axes[0]["Kp"] = YODA::Axis<double>({0.03, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19});
      _axes[0]["p"] = YODA::Axis<double>({0.06, 0.14, 0.155, 0.185, 0.215, 0.245, 0.275});
      _axes[0]["pi0"] = YODA::Axis<double>({0.1, 0.2, 0.3, 0.4, 0.5});
      _axes[0]["K0"] = YODA::Axis<double>({0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                           0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9});
      _axes[0]["lam"] = YODA::Axis<double>({0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                                            0.5, 0.65, 0.8, 0.95});
      _axes[0]["xi"] = YODA::Axis<double>({0.2, 0.3, 0.4, 0.5, 0.6, 0.7});
      _axes[0]["rho"] = YODA::Axis<double>({0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});
      _axes[0]["Kstarp"] = YODA::Axis<double>({0.06, 0.12, 0.24, 0.36, 0.48, 0.6});
      _axes[0]["Kstar0"] = YODA::Axis<double>({0.0, 0.06, 0.12, 0.24, 0.36, 0.48});
      _axes[0]["phi"] = YODA::Axis<double>({0.195, 0.385, 0.575, 0.945});

      _axes[1]["pip"] = YODA::Axis<double>({0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.58, 0.68, 0.88});
      _axes[1]["Kp"] = YODA::Axis<double>({0.02, 0.1, 0.11, 0.13, 0.15, 0.17, 0.19});
      _axes[1]["p"] = _axes[0]["p"];
      _axes[1]["pi0"] = _axes[0]["pi0"];
      _axes[1]["K0"] = _axes[0]["K0"];
      _axes[1]["lam"] = _axes[0]["lam"];
      _axes[1]["xi"] = _axes[0]["xi"];
      _axes[1]["rho"] = YODA::Axis<double>({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7});
      _axes[1]["Kstarp"] = _axes[0]["Kstarp"];
      _axes[1]["Kstar0"] = YODA::Axis<double>({0.06, 0.12, 0.24, 0.36, 0.48});
      _axes[1]["phi"] = YODA::Axis<double>({0.28, 0.36, 0.7, 1.0});

    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.abspid();
        if (id == PID::PIPLUS || id==PID::KPLUS   || id==PID::PROTON ||
            id==PID::PI0      || id==PID::K0S     || id==PID::K0L    ||
            id==PID::LAMBDA   || id==PID::XIMINUS || id==PID::RHO0   ||
            id==323 || id==313 || id==225 || id==PID::PHI) {
          unstable.push_back(p);
        }
        if(!p.children().empty())
          findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges[0].empty()) {
        for (const auto& item : _cont) {
          _edges[0][item.first] = item.second->xEdges();
          _edges[1][item.first] = _ups1[item.first]->xEdges();
        }
      }
      // Find the upsilons
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      // continuum
      if (upsilons.empty()) {
        _weightSum_cont->fill();
        const FinalState& fs = apply<FinalState>(event, "FS");
        // FS particles
        for (const Particle& p : fs.particles()) {
          int id = p.abspid();
          double xp = 2.*p.p3().mod()/sqrtS();
          if(id==PID::PIPLUS) {
            discfill("pip", xp, 0);
            _mult[1][0]->fill();
          }
          else if(id==PID::KPLUS) {
            discfill("Kp", xp, 0);
            _mult[1][1]->fill();
          }
          else if(id==PID::PROTON) {
            discfill("p", xp, 0);
            _mult[1][2]->fill();
          }
        }
        // Unstable particles
        for (const Particle& p : ufs.particles()) {
          int id = p.abspid();
          double xp = 2.*p.p3().mod()/sqrtS();
          if(id==PID::PI0) {
            discfill("pi0", xp, 0);
            _mult[1][3]->fill();
          }
          else if(id==PID::K0S || id==PID::K0L) {
            discfill("K0", xp, 0);
            _mult[1][4]->fill();
          }
          else if(id==PID::LAMBDA) {
            discfill("lam", xp, 0);
            _mult[1][5]->fill();
          }
          else if(id==PID::XIMINUS) {
            discfill("xi", xp, 0);
            _mult[1][6]->fill();
          }
          else if(id==PID::RHO0) {
            discfill("rho", xp, 0);
            _mult[1][7]->fill();
          }
          else if(id==323) {
            discfill("Kstarp", xp, 0);
            _mult[1][8]->fill();
          }
          else if(id==313) {
            discfill("Kstar0", xp, 0);
            _mult[1][9]->fill();
          }
          else if(id==PID::PHI) {
            discfill("phi", xp, 0);
            _mult[1][10]->fill();
          }
          else if(id==225) {
            _mult[1][11]->fill();
          }
        }
      }
      else {
        for (const Particle& ups : upsilons) {
          _weightSum_Ups1->fill();
          Particles unstable;
          LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          // Find the decay products we want
          findDecayProducts(ups,unstable);
          for (const Particle& p : unstable)  {
            int id = p.abspid();
            double xp = 2.*boost.transform(p.momentum()).p3().mod()/ups.mass();
            if(id==PID::PIPLUS) {
              discfill("pip", xp, 1);
              _mult[0][0]->fill();
            }
            else if(id==PID::KPLUS) {
              discfill("Kp", xp, 1);
              _mult[0][1]->fill();
            }
            else if(id==PID::PROTON) {
              discfill("p", xp, 1);
              _mult[0][2]->fill();
            }
            else if(id==PID::PI0) {
              discfill("pi0", xp, 1);
              _mult[0][3]->fill();
            }
            else if(id==PID::K0S || id==PID::K0L) {
              discfill("K0", xp, 1);
              _mult[0][4]->fill();
            }
            else if(id==PID::LAMBDA) {
              discfill("lam", xp, 1);
              _mult[0][5]->fill();
            }
            else if(id==PID::XIMINUS) {
              discfill("xi", xp, 1);
              _mult[0][6]->fill();
            }
            else if(id==PID::RHO0) {
              discfill("rho", xp, 1);
              _mult[0][7]->fill();
            }
            else if(id==323) {
              discfill("Kstarp", xp, 1);
              _mult[0][8]->fill();
            }
            else if(id==313) {
              discfill("Kstar0", xp, 1);
              _mult[0][9]->fill();
            }
            else if(id==PID::PHI) {
              discfill("phi", xp, 1);
              _mult[0][10]->fill();
            }
            else if(id==225) {
              _mult[0][11]->fill();
            }
          }
        }
      }
    }

    void discfill(const string& name, const double value, const size_t k) {
      string edge = "OTHER";
      const size_t idx = _axes[k][name].index(value);
      if (idx && idx <= _edges[k][name].size())  edge = _edges[k][name][idx-1];
      (k? _ups1 : _cont)[name]->fill(edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // multiplicities
      const vector<CounterPtr> scales = {_weightSum_Ups1, _weightSum_cont};
      for (size_t ix=0; ix<12; ++ix) {
        BinnedEstimatePtr<string> est;
        book(est, ix+12, 1, 1);
        for (size_t iy=0; iy<2; ++iy) {
          if (scales[iy]->val() > 0.) {
            unsigned int iz = iy==0 ? 2 : 1;
            scale(_mult[iy][ix], 1./ *scales[iy]);
            est->bin(iz).set(_mult[iy][ix]->val(), _mult[iy][ix]->err());
          }
        }
      }
      // spectra
      if (_weightSum_cont->val() > 0.) {
        scale(_cont, 1. / *_weightSum_cont);
        for( auto & hist : _cont) {
          for(auto & b: hist.second->bins()) {
            const size_t idx = b.index();
            b.scaleW(1./_axes[0][hist.first].width(idx));
          }
        }
      }
      if (_weightSum_Ups1->val() > 0.) {
        scale(_ups1, 1. / *_weightSum_Ups1);
        for( auto & hist : _ups1) {
          for(auto & b: hist.second->bins()) {
            const size_t idx = b.index();
            b.scaleW(1./_axes[1][hist.first].width(idx));
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,BinnedHistoPtr<string>> _cont, _ups1;
    map<string, YODA::Axis<double>> _axes[2];
    map<string, vector<string>> _edges[2];
    CounterPtr _weightSum_cont, _weightSum_Ups1;
    CounterPtr _mult[2][12];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1985_I205668);

}
