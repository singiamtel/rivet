// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Mass distribution in $e^+e^-\to e^+e^-D\bar{D}$ via $\gamma\gamma\to D\bar{D}$
  class BELLE_2006_I700451 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I700451);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(Cuts::abspid==411 || Cuts::abspid==421), "UFS");
      // histos
      for (unsigned int ix=0;ix<3;++ix) {
        book(_h_mass[ix],1,1,1+ix);
        if (ix<2) book(_h_mass_cut[ix],2,1,1+ix);
      }
      book(_h_angle,3,1,1);
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int &ncount) {
      for (const Particle& child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    bool findScattered(Particle beam, double& q2) {
      bool found = false;
      Particle scat = beam;
      while (!scat.children().empty()) {
        found = false;
        for (const Particle & p : scat.children()) {
          if (p.pid()==scat.pid()) {
            scat=p;
            found=true;
            break;
          }
        }
        if (!found) break;
      }
      if (!found) return false;
      q2 = -(beam.mom() - scat.mom()).mass2();
      return true;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find scattered leptons and calc Q2
      const Beam& beams = apply<Beam>(event, "Beams");
      double q12 = -1, q22 = -1;
      if (!findScattered(beams.beams().first,  q12)) vetoEvent;
      if (!findScattered(beams.beams().second, q22)) vetoEvent;
      // check the final state
      const FinalState & fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find the meson
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p1 : ufs.particles()) {
        if (p1.children().empty() || p1.pid()<0) continue;
        bool matched=false;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        for (const Particle& p2 : ufs.particles(Cuts::pid==-p1.pid())) {
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          matched = true;
          for (const auto& val : nRes2) {
            if (abs(val.first)==11) {
              if (val.second!=1) {
                matched = false;
                break;
              }
            }
            else if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            FourMomentum pDD = p1.mom()+p2.mom();
            double mDD = pDD.mass();
            if (p1.abspid()==421) _h_mass[0]->fill(mDD);
            else                  _h_mass[1]->fill(mDD);
            _h_mass[2]->fill(mDD);
            LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(pDD.betaVec());
            FourMomentum pD = boost.transform(p1.mom());
            double cTheta = abs(pD.z()/pD.p3().mod());
            if (abs(cTheta)<0.5) _h_mass_cut[0]->fill(mDD);
            else                 _h_mass_cut[1]->fill(mDD);
            if (mDD>3.91 && mDD<3.95) _h_angle->fill(cTheta);
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
	    normalize(_h_mass, 1.0, false);
	    normalize(_h_mass_cut, 1.0, false);
      normalize(_h_angle, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr  _h_mass[3],_h_mass_cut[2],_h_angle;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2006_I700451);

}
