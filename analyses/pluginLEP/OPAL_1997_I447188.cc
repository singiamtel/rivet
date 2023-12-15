// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Lambda polarization at LEP1
  class OPAL_1997_I447188 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1997_I447188);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const ChargedFinalState cfs;
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(UnstableParticles(), "UFS");

      // Book the histograms
      const vector<double> edges{0.027, 0.05, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 1.0};
      book(_h_ctheta, edges); size_t iy = 1;
      for (auto& b : _h_ctheta->bins()) {
        if (b.index() == 2 || b.index() == 5 || b.index() == 7) {
          book(b, 4, 1, iy++);
        }
        else {
          book(b , "/TMP/ctheta_"+std::to_string(b.index()-1), 20 , -1.0, 1.0);
        }
      }
      book(_h_ctheta_large, 4, 1, 4);

      book(_h_cphi, {0.3, 0.6, 0.9, 1.2, 1.5});
      //book(_h_cphi->bin(1), "/TMP/cphiP_0", 10, 0.0, 1.0);
      book(_h_cphi->bin(1), 5, 1, 1);
      book(_h_cphi->bin(2), 5, 1, 2);
      book(_h_cphi->bin(3), "/TMP/cphiP_3", 10, 0.0, 1.0);
      book(_h_cphi->bin(4), "/TMP/cphiP_4", 10, 0.0, 1.0);

      book(_h_cphi_low, 5, 1, 4);
      book(_h_cphi_mid, "/TMP/cphiP_mid", 10, 0.0, 1.0);
      book(_h_cphi_high, 5, 1, 3);

      book(_h_plus_lam, edges);
      book(_h_minus_lam, edges);
      for (size_t ix=0; ix<_h_plus_lam->numBins(); ++ix) {
        book(_h_plus_lam->bin(ix+1),  "/TMP/lamP_"+std::to_string(ix), 20, -1.0, 1.0);
        book(_h_minus_lam->bin(ix+1), "/TMP/lamM_"+std::to_string(ix), 20, -1.0, 1.0);
      }
      book(_h_plus_lam_large1, "/TMP/lamP_large_1",20,-1.,1.);
      book(_h_plus_lam_large2, "/TMP/lamP_large_2",20,-1.,1.);
      book(_h_minus_lam_large1, "/TMP/lamM_large_1",20,-1.,1.);
      book(_h_minus_lam_large2, "/TMP/lamM_large_2",20,-1.,1.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      Vector3 beamAxis;
      if (beams.first.pid()==-11) {
        beamAxis = beams.first .momentum().p3().unit();
      }
      else {
        beamAxis = beams.second.momentum().p3().unit();
      }

      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // thrust, to define an axis
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle & lambda : ufs.particles(Cuts::abspid==3122)) {
        double xE = lambda.momentum().t()/meanBeamMom;
        int sign = lambda.pid()/3122;
        Vector3 axis1 = lambda.momentum().p3().unit();
        // assymetry
        double cLam = axis1.dot(beamAxis);
        if (sign>0)
          _h_plus_lam->fill(xE,cLam);
        else {
          _h_minus_lam->fill(xE,cLam);
        }
        if(xE>0.15) {
          if (sign>0)
            _h_plus_lam_large1 ->fill(cLam);
          else {
            _h_minus_lam_large1->fill(cLam);
          }
        }
        if (xE>0.3) {
          if (sign>0)
            _h_plus_lam_large2 ->fill(cLam);
          else {
            _h_minus_lam_large2->fill(cLam);
          }
        }
        if(lambda.children().size()!=2) continue;
        // look at the decay products
        Particle proton,pion;
        if (lambda.children()[0].pid()==sign*2212 && lambda.children()[1].pid()==-sign*211) {
          proton = lambda.children()[0];
          pion   = lambda.children()[1];
        }
        else if (lambda.children()[1].pid()==sign*2212 && lambda.children()[0].pid()==-sign*211) {
          proton = lambda.children()[1];
          pion   = lambda.children()[0];
        }
        else {
          continue;
        }
        // boost to the lambda rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(lambda.momentum().betaVec());
        FourMomentum pproton = boost.transform(proton.momentum());
        // longitudinal polarization
        double ctheta = axis1.dot(pproton.p3().unit());
        _h_ctheta->fill(xE,ctheta);
        if(xE>=0.3)  _h_ctheta_large->fill(ctheta);
        // transverse polarization
        Vector3 axis2;
        if (lambda.momentum().p3().dot(thrust.thrustAxis())>=0.) {
          axis2 = thrust.thrustAxis();
        }
        else {
          axis2 =-thrust.thrustAxis();
        }
        Vector3 axis3 = axis2.cross(axis1).unit();
        double pT = sqrt(sqr(thrust.thrustMajorAxis().dot(lambda.momentum().p3()))+
                    sqr(thrust.thrustMinorAxis().dot(lambda.momentum().p3())));
        double cPhi = axis3.dot(pproton.p3().unit());
        _h_cphi->fill(pT,cPhi);
        if(pT>0.3) _h_cphi_low->fill(cPhi);
        if(pT>0.6) _h_cphi_mid->fill(cPhi);
        if(pT>1.5) _h_cphi_high->fill(cPhi);
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins()) {
        double Oi = bin.sumW();
        if(Oi==0.) continue;
        double ai = 0.5*(bin.xMax()-bin.xMin());
        double bi = 0.5*ai*(bin.xMax()+bin.xMin());
        double Ei = bin.errW();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    pair<double,double> calcAsymmetry(Estimate1DPtr hist,unsigned int mode) {
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins()) {
        double Oi = bin.val();
        if(Oi==0.) continue;
        double bi;
        if(mode==0)
          bi = 0.25*(bin.xMax()-bin.xMin())*(bin.xMax()+bin.xMin());
        else
          bi = 4.*(bin.xMax()+bin.xMin())/(3.+sqr(bin.xMax())+bin.xMax()*bin.xMin()+sqr(bin.xMin()));
        double Ei = bin.errAvg();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*Oi;
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // longitudinal polarization
      double aLam = 0.642;
      Estimate1DPtr h_long;
      book(h_long, 1, 1, 1);
      for (auto& hist : _h_ctheta->bins()) {
        normalize(hist);
        pair<double,double> alpha = calcAlpha(hist);
        alpha.first  /=aLam;
        alpha.second /=aLam;
        h_long->bin(hist.index()).set(100.*alpha.first, 100.*alpha.second);
      }
      normalize(_h_ctheta_large);
      pair<double,double> alpha = calcAlpha(_h_ctheta_large);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_long_l;
      book(h_long_l,1,1,2);
      h_long_l->bin(1).set(100.*alpha.first, 100.*alpha.second);
      // transverse polarization
      Estimate1DPtr h_trans;
      book(h_trans,2,1,1);
      for (auto& hist : _h_cphi->bins()) {
        normalize(hist);
        pair<double,double> alpha = calcAlpha(hist);
        alpha.first  /=aLam;
        alpha.second /=aLam;
        h_trans->bin(hist.index()).set(100.*alpha.first, 100.*alpha.second);
      }
      normalize(_h_cphi_low);
      alpha = calcAlpha(_h_cphi_low);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_trans_low;
      book(h_trans_low,2,1,2);
      h_trans_low->bin(1).set(alpha.first, 100.*alpha.second);
      normalize(_h_cphi_mid);
      alpha = calcAlpha(_h_cphi_mid);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_trans_mid;
      book(h_trans_mid,2,1,3);
      h_trans_mid->bin(1).set(alpha.first, 100.*alpha.second);
      normalize(_h_cphi_high);
      alpha = calcAlpha(_h_cphi_high);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_trans_high;
      book(h_trans_high,2,1,4);
      h_trans_high->bin(1).set(alpha.first, 100.*alpha.second);
      // asyymetry
      Estimate1DPtr h_asym;
      book(h_asym, 3, 1, 1);
      for (size_t ix=0; ix < _h_plus_lam->numBins(); ++ix) {
       	normalize(_h_plus_lam->bin(ix+1));
       	normalize(_h_minus_lam->bin(ix+1));
       	Estimate1DPtr eTemp;
        book(eTemp, "/TMP/a_lam_" + to_string(ix));
       	asymm(_h_plus_lam->bin(ix+1), _h_minus_lam->bin(ix+1), eTemp);
       	pair<double,double> alpha = calcAsymmetry(eTemp, 1);
       	h_asym->bin(ix+1).set(-alpha.first, alpha.second);
      }
      normalize(_h_plus_lam_large1 );
      normalize(_h_minus_lam_large1);
      Estimate1DPtr eTemp;
      book(eTemp, "/TMP/a_lam_large1");
      asymm(_h_plus_lam_large1,_h_minus_lam_large1,eTemp);
      alpha = calcAsymmetry(eTemp,1);
      book(h_asym,3,1,2);
      h_asym->bin(1).set(-alpha.first, alpha.second);
      normalize(_h_plus_lam_large2 );
      normalize(_h_minus_lam_large2);
      book(eTemp, "/TMP/a_lam_large2");
      asymm(_h_plus_lam_large2,_h_minus_lam_large2, eTemp);
      alpha = calcAsymmetry(eTemp,1);
      book(h_asym,3,1,3);
      h_asym->bin(1).set(-alpha.first, alpha.second);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_ctheta,_h_cphi,_h_plus_lam,_h_minus_lam;
    Histo1DPtr _h_ctheta_large;
    Histo1DPtr _h_cphi_low, _h_cphi_mid, _h_cphi_high;
    Histo1DPtr _h_plus_lam_large1,_h_plus_lam_large2;
    Histo1DPtr _h_minus_lam_large1,_h_minus_lam_large2;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(OPAL_1997_I447188);


}
