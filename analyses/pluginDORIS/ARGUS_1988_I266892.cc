// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief baryon correlations
  class ARGUS_1988_I266892 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1988_I266892);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs                 , "FS"    );
      declare(UnstableParticles(), "UFS"   );
      declare(Thrust(fs)         , "Thrust");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_nBB2[ix], "TMP/nBB2_"+toString(ix+1),
             refData<YODA::BinnedEstimate<string>>(4+2*ix, 1, 1));
        book(_n[ix],    "TMP/n_"+toString(ix+1),
             refData<YODA::BinnedEstimate<string>>(4+2*ix, 1, 1));
        book(_nBB [ix],3+2*ix,1,1);
        book(_c[ix],"TMP/c_"+toString(ix+1));
        if (ix==2) continue;
        book(_h_pt[ix],1,1,1+ix);
      }
      book(_h_pt[2],2,1,1);
    }

    /// Recursively walk the decay tree to find the stable decay products of @a p
    void findDecayProducts(const Particle& mother, Particles& final,
                           Particles& protons, Particles& lambda, Particles& xi,
                           Particles& l1520) {
      for (const Particle & p: mother.children()) {
        if      (p.abspid()==3122  ) lambda.push_back(p);
        else if (p.abspid()==3312  ) xi    .push_back(p);
        else if (p.abspid()==102134) l1520 .push_back(p);
        if (!p.children().empty()) {
          findDecayProducts(p, final, protons, lambda, xi, l1520);
        }
        else {
          final.push_back(p);
          if (p.abspid()==2212) protons.push_back(p);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_edges.empty()) _edges = _nBB2[1]->xEdges();
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      if (upsilons.empty()) {
        _c[1]->fill();
        _c[2]->fill();
        Thrust thrust = apply<Thrust>(event, "Thrust");
        Vector3 axis = thrust.thrustAxis();
        Particles protons = apply<FinalState>(event, "FS").particles(Cuts::abspid==2212);
        for (unsigned int ix=0;ix<protons.size();++ix) {
          double d1   = protons[ix].mom().p3().dot(axis);
          Vector3 pt1 = protons[ix].mom().p3()-d1*axis;
          for (unsigned int iy=ix+1;iy<protons.size();++iy) {
            if (protons[ix].pid()*protons[iy].pid()>0) continue;
            double d2 = protons[iy].mom().p3().dot(axis);
            Vector3 pt2 = protons[iy].mom().p3()-d2*axis;
            double phi = pt1.angle(pt2);
            if (phi<0.) phi +=M_PI;
            phi *=180./M_PI;
            if(phi==180.) phi=phi-1e-10;
            if (d1*d2>0) _h_pt[0]->fill(phi);
            else         _h_pt[2]->fill(phi);
          }
        }
        Particles lambda = ufs.particles(Cuts::abspid==3122);
        for (unsigned int ix=0; ix<lambda.size(); ++ix) {
          _n[0]->fill(_edges[1]);
          for (unsigned int iy=ix+1; iy<lambda.size(); ++iy) {
            if (lambda[ix].pid()*lambda[iy].pid()<0) {
              _nBB [0]->fill(_edges[1]);
              _nBB2[0]->fill(_edges[1]);
            }
          }
        }
        Particles xi = ufs.particles(Cuts::abspid==3312);
        for (unsigned int ix=0;ix<xi.size();++ix) {
          _n[1]->fill(_edges[1]);
          _n[1]->fill(_edges[2]);
          for (unsigned int iy=0;iy<lambda.size();++iy) {
            if (xi[ix].pid()*lambda[iy].pid()<0) {
              _nBB [1]->fill(_edges[1]);
              _nBB2[1]->fill(_edges[1]);
              _nBB2[1]->fill(_edges[2]);
            }
          }
        }
        Particles l1520 = ufs.particles(Cuts::abspid==102134);
        for (unsigned int ix=0; ix<l1520.size(); ++ix) {
          _n[2]->fill(_edges[1]);
          _n[2]->fill(_edges[2]);
          for (unsigned int iy=0; iy<lambda.size(); ++iy) {
            if (l1520[ix].pid()*lambda[iy].pid()<0) {
              _nBB [2]->fill(_edges[1]);
              _nBB2[2]->fill(_edges[1]);
              _nBB2[2]->fill(_edges[2]);
            }
          }
        }
      }
      else {
        for (const Particle& ups : upsilons) {
          _c[0]->fill();
          _c[2]->fill();
          LorentzTransform boost;
          if (ups.p3().mod() > 1*MeV) {
            boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
          }
          Particles final, protons, lambda, xi, l1520;
          findDecayProducts(ups, final, protons, lambda, xi, l1520);
          vector<FourMomentum> mom; mom.reserve(final.size());
          for (const Particle& p : final) {
            mom.push_back(boost.transform(p.mom()));
          }
          Thrust thrust;
          thrust.calc(mom);
          Vector3 axis = thrust.thrustAxis();
          for (unsigned int ix=0; ix<protons.size(); ++ix) {
            double d1   = protons[ix].mom().p3().dot(axis);
            Vector3 pt1 = protons[ix].mom().p3()-d1*axis;
            for (unsigned int iy=ix+1;iy<protons.size();++iy) {
              if (protons[ix].pid()*protons[iy].pid()>0) continue;
              double d2   = protons[iy].mom().p3().dot(axis);
              Vector3 pt2 = protons[iy].mom().p3()-d2*axis;
              double phi = pt1.angle(pt2);
              if (phi<0.) phi +=M_PI;
              phi *=180./M_PI;
              if(phi==180.) phi=phi-1e-10;
              if (d1*d2>0) _h_pt[1]->fill(phi);
            }
          }
          for (unsigned int ix=0; ix<lambda.size(); ++ix) {
            _n[0]->fill(_edges[0]);
            for (unsigned int iy=ix+1; iy<lambda.size(); ++iy) {
              if (lambda[ix].pid()*lambda[iy].pid()<0) {
                _nBB [0]->fill(_edges[0]);
                _nBB2[0]->fill(_edges[0]);
              }
            }
          }
          for (unsigned int ix=0; ix<xi.size(); ++ix) {
            _n[1]->fill(_edges[0]);
            _n[1]->fill(_edges[2]);
            for (unsigned int iy=0; iy<lambda.size(); ++iy) {
              if(xi[ix].pid()*lambda[iy].pid()<0) {
                _nBB [1]->fill(_edges[0]);
                _nBB2[1]->fill(_edges[0]);
                _nBB2[1]->fill(_edges[2]);
              }
            }
          }
          for (unsigned int ix=0; ix<l1520.size(); ++ix) {
            _n[2]->fill(_edges[0]);
            _n[2]->fill(_edges[2]);
            for (unsigned int iy=0; iy<lambda.size(); ++iy) {
              if (l1520[ix].pid()*lambda[iy].pid()<0) {
                _nBB [2]->fill(_edges[0]);
                _nBB2[2]->fill(_edges[0]);
                _nBB2[2]->fill(_edges[2]);
              }
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pt, 1.0, false);
      for (unsigned int ix=0; ix<3; ++ix) {
        for(unsigned int iy=0;iy<_nBB[ix]->numBins();++iy)
          _nBB[ix]->bin(iy+1).scaleW( 1./_c[iy]->val());
        BinnedEstimatePtr<string> tmp;
        book(tmp,4+2*ix, 1, 1);
        divide(_nBB2[ix], _n[ix],tmp);
        if(ix==0) scale(tmp,2.);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pt[3];
    BinnedHistoPtr<string> _nBB[3],_nBB2[3],_n[3];
    CounterPtr _c[3];
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1988_I266892);

}
