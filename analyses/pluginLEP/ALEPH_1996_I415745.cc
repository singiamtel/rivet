// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Lambda polarization at LEP1
  class ALEPH_1996_I415745 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_1996_I415745);


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

      // Book histograms
      book(_h_ctheta, {0.1, 0.15, 0.2, 0.3, 0.4, 1.});
      for (size_t i = 0; i < _h_ctheta->numBins(); ++i) {
        book(_h_ctheta->bin(i+1), "/TMP/ctheta_"+to_string(i), 20, -1.0, 1.0);
      }
      book(_h_ctheta_large,"/TMP/ctheta_large",20,-1.,1.);

      book(_h_plus_cphi, {0.3, 0.6, 0.9, 1.2, 1.5});
      book(_h_minus_cphi, {0.3, 0.6, 0.9, 1.2, 1.5});
      for (size_t i = 0; i < _h_plus_cphi->numBins(); ++i) {
        book(_h_plus_cphi->bin(i+1), "/TMP/cphiP_0_"+to_string(i), 10, 0., 1.);
        book(_h_minus_cphi->bin(i+1), "/TMP/cphiM_0_"+to_string(i), 10, 0., 1.);
      }
      book(_h_plus_cphi_low  , "/TMP/cphiP_low" ,10,0.,1.);
      book(_h_plus_cphi_mid  , "/TMP/cphiP_mid" ,10,0.,1.);
      book(_h_plus_cphi_high , "/TMP/cphiP_high",10,0.,1.);
      book(_h_minus_cphi_low , "/TMP/cphiM_low" ,10,0.,1.);
      book(_h_minus_cphi_mid , "/TMP/cphiM_mid" ,10,0.,1.);
      book(_h_minus_cphi_high, "/TMP/cphiM_high",10,0.,1.);

      book(_h_plus_lam, {0.1, 0.15, 0.2, 0.3, 0.4, 1.});
      book(_h_minus_lam, {0.1, 0.15, 0.2, 0.3, 0.4, 1.});
      for (size_t i = 0; i < _h_plus_lam->numBins(); ++i) {
        book(_h_plus_lam->bin(i+1), "/TMP/lamP_0_"+to_string(i), 20, -1., 1.);
        book(_h_minus_lam->bin(i+1), "/TMP/lamM_0_"+to_string(i), 20, -1., 1.);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
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
        double z = lambda.momentum().p3().mod()/meanBeamMom;
        int sign = lambda.pid()/3122;
        Vector3 axis1 = lambda.momentum().p3().unit();
        // assymetry
        double cLam = axis1.dot(beamAxis);
        if(sign>0)
          _h_plus_lam->fill(z,cLam);
        else
          _h_minus_lam->fill(z,cLam);
        if(lambda.children().size()!=2) continue;
        // look at the decay products
        Particle proton,pion;
        if(lambda.children()[0].pid()==sign*2212 &&
           lambda.children()[1].pid()==-sign*211) {
          proton = lambda.children()[0];
          pion   = lambda.children()[1];
        }
        else if(lambda.children()[1].pid()==sign*2212 &&
          lambda.children()[0].pid()==-sign*211) {
          proton = lambda.children()[1];
          pion   = lambda.children()[0];
        }
        else
          continue;
        // boost to the lambda rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(lambda.momentum().betaVec());
        FourMomentum pproton = boost.transform(proton.momentum());
        // longitudinal polarization
        double ctheta = axis1.dot(pproton.p3().unit());
        _h_ctheta->fill(z,ctheta);
        if(z>=0.3)  _h_ctheta_large->fill(ctheta);

        // transverse polarization
        if (z>0.15) {
          Vector3 axis2;
          if(lambda.momentum().p3().dot(thrust.thrustAxis())>=0.) {
            axis2 = thrust.thrustAxis();
          }
          else {
            axis2 =-thrust.thrustAxis();
          }
          Vector3 axis3 = axis2.cross(axis1).unit();

          double pT = sqrt(sqr(thrust.thrustMajorAxis().dot(lambda.momentum().p3()))+
               sqr(thrust.thrustMinorAxis().dot(lambda.momentum().p3())));
          double cPhi = axis3.dot(pproton.p3().unit());
          if(cPhi>0.) {
            _h_plus_cphi->fill(pT,cPhi);
            if(pT>0.3)
              _h_plus_cphi_low->fill(cPhi);
            if(pT>0.6)
              _h_plus_cphi_mid->fill(cPhi);
            if(pT>1.5)
              _h_plus_cphi_high->fill(cPhi);
          }
          else {
            _h_minus_cphi->fill(pT,abs(cPhi));
            if(pT>0.3)
              _h_minus_cphi_low->fill(abs(cPhi));
            if(pT>0.6)
              _h_minus_cphi_mid->fill(abs(cPhi));
            if(pT>1.5)
              _h_minus_cphi_high->fill(abs(cPhi));
          }
        }
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins() ) {
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

    pair<double,double> calcAsymmetry(Estimate1DPtr hist, unsigned int mode) {
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins() ) {
        double Oi = bin.val();
        if(Oi==0.) continue;
        double bi;
        if (mode==0)
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
      unsigned int ipoint=0;
      double aLam = 0.642;
      Estimate1DPtr h_long;
      book(h_long,1,1,1);
      for (auto& hist : _h_ctheta->bins()) {
        normalize(hist);
        pair<double,double> alpha = calcAlpha(hist);
        alpha.first  /=aLam;
        alpha.second /=aLam;
        h_long->bin(ipoint+1).set(alpha.first, alpha.second);
        ++ipoint;
      }
      normalize(_h_ctheta_large);
      pair<double,double> alpha = calcAlpha(_h_ctheta_large);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_long_l;
      book(h_long_l,1,2,1);
      h_long_l->bin(1).set(alpha.first, alpha.second);
      // transverse polarization
      Estimate1DPtr h_trans;
      book(h_trans,2,1,1);
      for (size_t ix=0; ix < _h_plus_cphi->numBins(); ++ix) {
        normalize(_h_plus_cphi->bin(ix));
        normalize(_h_minus_cphi->bin(ix));
        Estimate1DPtr sTemp;
        book(sTemp, "/TMP/a_cphi_"+to_string(ix),10,0.,1.);
        asymm(_h_plus_cphi->bin(ix+1), _h_minus_cphi->bin(ix+1), sTemp);
        pair<double,double> alpha = calcAsymmetry(sTemp,0);
        alpha.first  /=aLam;
        alpha.second /=aLam;
        h_trans->bin(ix+1).set(alpha.first, alpha.second);
      }
      Estimate1DPtr sLow;
      book(sLow,"/TMP/a_cphi_low",10,0.,1.);
      asymm(_h_plus_cphi_low, _h_minus_cphi_low, sLow);
      alpha = calcAsymmetry(sLow,0);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_trans_low;
      book(h_trans_low,2,3,1);
      h_trans_low->bin(1).set(alpha.first, alpha.second);

      Estimate1DPtr sMid;
      book(sMid,"/TMP/a_cphi_mid",10,0.,1.);
      asymm(_h_plus_cphi_mid,_h_minus_cphi_mid,sMid);
      alpha = calcAsymmetry(sMid,0);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_trans_mid;
      book(h_trans_mid,2,4,1);
      h_trans_mid->bin(1).set(alpha.first, alpha.second);

      Estimate1DPtr sHigh;
      book(sHigh,"/TMP/a_cphi_high",10,0.,1.);
      asymm(_h_plus_cphi_high,_h_minus_cphi_high,sHigh);
      alpha = calcAsymmetry(sHigh,0);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Estimate1DPtr h_trans_high;
      book(h_trans_high,2,2,1);
      h_trans_high->bin(1).set(alpha.first, alpha.second);

      // asyymetry
      Estimate1DPtr h_asym;
      book(h_asym,3,1,1);
      for (size_t ix=0; ix < _h_plus_lam->numBins(); ++ix) {
        normalize(_h_plus_lam->bin(ix+1));
        normalize(_h_minus_lam->bin(ix+1));
        Estimate1DPtr sTemp;
        book(sTemp, "/TMP/a_lam_"+to_string(ix), 20, -1., 1.);
        asymm(_h_plus_lam->bin(ix+1), _h_minus_lam->bin(ix+1), sTemp);
        pair<double,double> alpha = calcAsymmetry(sTemp,1);
        h_asym->bin(ix+1).set(alpha.first, alpha.second);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_ctheta,_h_plus_cphi,_h_minus_cphi,_h_plus_lam,_h_minus_lam;
    Histo1DPtr _h_ctheta_large;
    Histo1DPtr _h_minus_cphi_low, _h_minus_cphi_mid, _h_minus_cphi_high;
    Histo1DPtr _h_plus_cphi_low, _h_plus_cphi_mid, _h_plus_cphi_high;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALEPH_1996_I415745);


}
