// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Collins assymmetry
  class BESIII_2016_I1384778 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2016_I1384778);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(FinalState(Cuts::abspid==PID::PIPLUS), "FS");
      // book the histograms
      _h_L.resize(6);
      _h_U.resize(6);
      _h_C.resize(6);
      for(size_t ix=0; _h_L.size(); ++ix) {
        const string pre = "/TMP/h_z1z2_"+to_string(ix+1);
        book(_h_L[ix], pre+"_L", 20, 0., M_PI);
        book(_h_U[ix], pre+"_U", 20, 0., M_PI);
        book(_h_C[ix], pre+"_C", 20, 0., M_PI);
      }
      book(_h_pT_L, {0., 0.2, 0.3, 0.45, 0.8, 1.4});
      book(_h_pT_U, {0., 0.2, 0.3, 0.45, 0.8, 1.4});
      book(_h_pT_C, {0., 0.2, 0.3, 0.45, 0.8, 1.4});
      for (size_t ix=1; ix < _h_pT_L->numBins()+1; ++ix) {
        const string pre = "/TMP/h_pT_"+to_string(ix);
        book(_h_pT_L->bin(ix), pre+"_L", 20, 0.0, M_PI);
        book(_h_pT_U->bin(ix), pre+"_U", 20, 0.0, M_PI);
        book(_h_pT_C->bin(ix), pre+"_C", 20, 0.0, M_PI);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.pid()>0) {
        axis = beams.first .momentum().p3().unit();
      }
      else {
        axis = beams.second.momentum().p3().unit();
      }
      // loop over pions pair, using index to avoid double counting
      Particles pions = apply<FinalState>(event, "FS").particles();
      for (unsigned int i1=0; i1<pions.size(); ++i1) {
        const double x1=2.*pions[i1].momentum().t()/sqrtS();
        // cut on z1
        if (x1<0.2||x1>0.9) continue;
        // cos theta cut
        if (abs(cos(pions[i1].momentum().p3().polarAngle()))>0.93) continue;
        for(unsigned int i2=i1+1;i2<pions.size();++i2) {
          // cut on z2
          const double x2=2.*pions[i2].momentum().t()/sqrtS();
          if (x2<0.2||x2>0.9) continue;
          // cos theta cut
          if (abs(cos(pions[i2].momentum().p3().polarAngle()))>0.93) continue;
          // cut on opening angle (>120 degrees)
          if (pions[i1].momentum().p3().angle(pions[i2].momentum().p3())>2.*M_PI/3.) {
            continue;
          }
          Particle p1=pions[i1], p2=pions[i2];
          double z1(x1),z2(x2);
          // randomly order the particles
          if (rand()/static_cast<double>(RAND_MAX) < 0.5) {
            swap(p1,p2);
            swap(z1,z2);
          }
          // particle 2 defines the z axis
          Vector3 ez = p2.momentum().p3().unit();
          // beam and 2 define the plane (y is normal to plane)
          Vector3 ey = ez.cross(axis).unit();
          // x by cross product
          Vector3 ex = ey.cross(ez).unit();
          // phi
          double phi = ex.angle(p1.momentum().p3());
          // hists vs z1,z2
          unsigned int ibin=0;
          if (z1<=.3&&z2<=.3) {
            ibin=0;
          }
          else if (z1>0.5&&z2>0.5) {
            ibin=5;
          }
          else if (min(z1,z2)<=0.3) {
            if (max(z1,z2)>0.5) {
              ibin=2;
            }
            else {
              ibin=1;
            }
          }
          else {
            if (max(z1,z2)>0.5) {
              ibin=4;
            }
            else {
              ibin=3;
            }
          }
          _h_C[ibin]->fill(phi);
          if (p1.pid()==p2.pid()) {
            _h_L[ibin]->fill(phi);
          }
          else {
            _h_U[ibin]->fill(phi);
          }
          // hists vs pT
          double pPar2 = sqr(ez.dot(p1.momentum().p3()));
          double pPerp = sqrt(p1.momentum().p3().mod2()-pPar2);
          _h_pT_C->fill(pPerp,phi);
          if(p1.pid()==p2.pid()) {
            _h_pT_L->fill(pPerp,phi);
          }
          else {
            _h_pT_U->fill(pPerp,phi);
          }
        }
      }
    }

    pair<double,double> calcAsymmetry(Estimate1DPtr hist) {
      double sum1(0.),sum2(0.);
      for (const auto& b : hist->bins()) {
        double Oi = b.val();
        if(Oi==0. || std::isnan(Oi) ) continue;
        double ai = 1.;
        double bi = 0.5*(sin(2.*b.xMax())-sin(2.*b.xMin()))/b.xWidth();
        double Ei = b.errAvg();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1, sqrt(1./sum1));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // ratios
      Estimate1DPtr _h_z_UL,_h_z_UC;
      book(_h_z_UL,1,1,5);
      book(_h_z_UC,1,1,6);
      for (size_t ix=0; ix < _h_L.size(); ++ix) {
        normalize(_h_L[ix]);
        normalize(_h_U[ix]);
        normalize(_h_C[ix]);
        const string pre = "/TMP/R_z1z2_"+to_string(ix+1);
        Estimate1DPtr R1;
        book(R1, pre+"_UL");
        divide(_h_U[ix], _h_L[ix], R1);
        Estimate1DPtr R2;
        book(R2, pre+"_UC");
        divide(_h_U[ix], _h_C[ix], R2);
        pair<double,double> asym1 = calcAsymmetry(R1);
        _h_z_UL->bin(ix+1).set(asym1.first, asym1.second);
        pair<double,double> asym2 = calcAsymmetry(R2);
        _h_z_UC->bin(ix+1).set(asym2.first, asym2.second);
      }
      Estimate1DPtr _h_pT_UL,_h_pT_UC;
      book(_h_pT_UL,2,1,4);
      book(_h_pT_UC,2,1,5);
      for (size_t ix=1; ix < _h_pT_L->numBins()+1; ++ix) {
        normalize(_h_pT_L->bin(ix));
        normalize(_h_pT_U->bin(ix));
        normalize(_h_pT_C->bin(ix));
        const string pre = "/TMP/R_pT_"+to_string(ix);
        Estimate1DPtr R1;
        book(R1, pre+"_UL");
        divide(_h_pT_U->bin(ix), _h_pT_L->bin(ix), R1);
        Estimate1DPtr R2;
        book(R2, pre+"_UC");
        divide(_h_pT_U->bin(ix), _h_pT_C->bin(ix), R2);
        pair<double,double> asym1 = calcAsymmetry(R1);
        _h_pT_UL->bin(ix).set(asym1.first, asym1.second);
        pair<double,double> asym2 = calcAsymmetry(R2);
        _h_pT_UC->bin(ix).set(asym2.first, asym2.second);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    vector<Histo1DPtr> _h_L,_h_U,_h_C;
    Histo1DGroupPtr _h_pT_L,_h_pT_U,_h_pT_C;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2016_I1384778);

}
