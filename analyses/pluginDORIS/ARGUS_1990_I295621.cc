// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief baryon correlations
  class ARGUS_1990_I295621 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1990_I295621);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(Cuts::pid==2212), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      for (unsigned int ix=0;ix<2;++ix) {
        book(_c[ix],"TMP/c_"+toString(ix+1));
        book(_h_angle[ix], 9, 1, 1+ix);
      }
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_n_cont[ix], 1+ix, 1, 1);
        book(_n_Ups [ix], 4+ix, 1, 1);
      }
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(const Particle& mother, Particles& protons, Particles& lam) {
      for (const Particle& p: mother.children()) {
        const int id = abs(p.pid());
        if (id == 2212) {
          protons.push_back(p);
        }
        else if (id == 3122) {
          lam.push_back(p);
        }
        if (!p.children().empty()) {
          findDecayProducts(p, protons, lam);
        }
      }
    }

    void fillHistos(unsigned int imode, Particles& protons, Particles& lam, LorentzTransform& cms_boost) {
      Particles protons2;
      for (const Particle& p : protons) {
        const double modp = cms_boost.transform(p).p3().mod();
        if (modp<0.4 || modp>1.2) continue;
        if (imode==0) _n_cont[0]->fill("< 9.46"s);
        else          _n_Ups [0]->fill("9.460 + 10.023"s);
        protons2.push_back(p);
      }
      for (unsigned int ix=0; ix<protons2.size(); ++ix) {
        Vector3 axis1 = cms_boost(protons2[ix].mom()).p3().unit();
        if (protons2[ix].parents()[0].abspid()==3122) continue;
        for (unsigned int iy=ix+1;iy<protons2.size();++iy) {
          if (protons2[ix].pid()*protons2[iy].pid()<0) continue;
          if (protons2[iy].parents()[0].abspid()==3122) continue;
          Vector3 axis2 = cms_boost(protons2[iy].mom()).p3().unit();
          double cTheta = axis1.dot(axis2);
          if (imode==0) _n_cont[1]->fill("< 9.46"s);
          else         _n_Ups [1]->fill("9.460 + 10.023"s);
          _h_angle[imode]->fill(cTheta);
        }
      }
      for (unsigned int ix=0; ix<lam.size(); ++ix) {
        for (unsigned int iy=ix+1; iy<lam.size(); ++iy) {
          if (lam[ix].pid()*lam[iy].pid()<0) continue;
          if (imode==0) _n_cont[2]->fill("< 9.46"s);
          else          _n_Ups [2]->fill("9.460 + 10.023"s);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the upsilons
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 || Cuts::pid==100553);
      // continuum
      if (upsilons.empty()) {
        _c[0]->fill();
        Particles protons = apply<FinalState>(event,"FS").particles();
        Particles lam = ufs.particles(Cuts::abspid==3122);
        LorentzTransform cms_boost;
        fillHistos(0,protons,lam,cms_boost);
      }
      // found an upsilon
      else {
        for (const Particle& ups : upsilons) {
          _c[1]->fill();
          Particles protons,lam;
          findDecayProducts(ups, protons, lam);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 0.001) {
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
          }
          fillHistos(1,protons,lam,cms_boost);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_n_cont, 1.0/ *_c[0]);
      scale(_n_Ups,  1.0/ *_c[1]);
      scale(_h_angle[0], 1e4/ *_c[0]);
      scale(_h_angle[1], 1e3/ *_c[1]);
      BinnedEstimatePtr<string> tmp;
      for (unsigned int ix=0;ix<2;++ix) {
        book(tmp, 7+ix, 1, 1);
        double val = _n_Ups[ix]->bin(1).sumW()/_n_cont[ix]->bin(1).sumW();
        double err = val*sqrt( sqr(_n_Ups [ix]->bin(1).errW()/ _n_Ups[ix]->bin(1).sumW()) +
                               sqr(_n_cont[ix]->bin(1).errW()/_n_cont[ix]->bin(1).sumW()));
        tmp->bin(1).set(val,err);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _n_cont[3],_n_Ups[3];
    Histo1DPtr _h_angle[2];
    CounterPtr _c[2];

    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1990_I295621);

}
