// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/Psi and Psi(2S) spectra at the Upsilon(4S)
  class CLEOII_2002_I606309 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2002_I606309);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      //histograms
      book(_weightSum,"/TMP/weightSum");
      book(_h_Jpsi            ,3,1,1);
      book(_h_Psi_prime       ,3,1,2);
      book(_h_cTheta_Jpsi     ,4,1,1);
      book(_h_cTheta_Psi_prime,4,1,2);
      book(_h2_cTheta_Jpsi, {0.,0.8,1.4,2.});
      for (auto& b : _h2_cTheta_Jpsi->bins()) {
        book(b, "/TMP/ctheta_"+to_string(b.index()-1), 20, -1.0, 1.0);
      }
    }

    void findDecayProducts(Particle parent, Particles& charm) {
      for (const Particle & p :parent.children()) {
        if (p.pid()==443) {
          charm.push_back(p);
          continue;
        }
        else if (p.pid()==100443) {
          charm.push_back(p);
        }
        if(!p.children().empty()) {
          findDecayProducts(p,charm);
        }
      }
    }

    void findLeptons(const Particle & mother, unsigned int& nstable, Particles& lp, Particles& lm) {
      for (const Particle& p : mother.children()) {
        int id = p.pid();
      	if ( id == 11 || id == 13 ) {
          lm.push_back(p);
          ++nstable;
        }
       	else if (id == -11 || id==-13) {
       	  lp.push_back(p);
       	  ++nstable;
       	}
        else if ( !p.children().empty() ) {
          findLeptons(p,nstable,lp,lm);
        }
        else {
          ++nstable;
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
        _weightSum->fill();
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        Particles charm;
        findDecayProducts(p,charm);
        for (const Particle& pp : charm) {
          FourMomentum pcm = boost.transform(pp.momentum());
          unsigned int nstable = 0;
          Particles lp, lm;
          findLeptons(pp,nstable,lp,lm);
          double cTheta(0.);
          bool foundLeptons(false);
          if (nstable==2&&lp.size()==1&&lm.size()==1) {
            foundLeptons=true;
            FourMomentum pl = boost.transform(lm[0].momentum());
            const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pcm.betaVec());
            pl = boost2.transform(pl);
            cTheta = pl.p3().unit().dot(pcm.p3().unit());
          }
          if(pp.pid()==443) {
            _h_Jpsi->fill(pcm.p3().mod());
            if(foundLeptons) {
              _h_cTheta_Jpsi->fill(cTheta);
              _h2_cTheta_Jpsi->fill(pcm.p3().mod(), cTheta);
            }
          }
          else {
            _h_Psi_prime->fill(pcm.p3().mod());
            if(foundLeptons) {
              _h_cTheta_Psi_prime->fill(cTheta);
            }
          }
        }
      }
    }

    pair<double,pair<double,double> > calcAlpha(const Histo1DPtr& hist) {
      if(hist->numEntries()==0.) return make_pair(0.,make_pair(0.,0.));
      double d = 3./(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double c = 3.*(hist->xMax()-hist->xMin())/(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double sum1(0.),sum2(0.),sum3(0.),sum4(0.),sum5(0.);
      for (const auto& bin : hist->bins() ) {
       	double Oi = bin.sumW();
        if(Oi==0.) continue;
        double a =  d*(bin.xMax() - bin.xMin());
        double b = d/3.*(pow(bin.xMax(),3) - pow(bin.xMin(),3));
       	double Ei = bin.errW();
        sum1 +=   a*Oi/sqr(Ei);
        sum2 +=   b*Oi/sqr(Ei);
        sum3 += sqr(a)/sqr(Ei);
        sum4 += sqr(b)/sqr(Ei);
        sum5 +=    a*b/sqr(Ei);
      }
      // calculate alpha
      double alpha = (-c*sum1 + sqr(c)*sum2 + sum3 - c*sum5)/(sum1 - c*sum2 + c*sum4 - sum5);
      // and error
      double cc = -pow((sum3 + sqr(c)*sum4 - 2*c*sum5),3);
      double bb = -2*sqr(sum3 + sqr(c)*sum4 - 2*c*sum5)*(sum1 - c*sum2 + c*sum4 - sum5);
      double aa =  sqr(sum1 - c*sum2 + c*sum4 - sum5)*(-sum3 - sqr(c)*sum4 + sqr(sum1 - c*sum2 + c*sum4 - sum5) + 2*c*sum5);
      double dis = sqr(bb)-4.*aa*cc;
      if (dis>0.) {
        dis = sqrt(dis);
        return make_pair(alpha,make_pair(0.5*(-bb+dis)/aa,-0.5*(-bb-dis)/aa));
      }
      else {
        return make_pair(alpha,make_pair(0.,0.));
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      if(_weightSum->val()==0.) return;
      scale(_h_Jpsi     , 50./ *_weightSum);
      scale(_h_Psi_prime, 50./ *_weightSum);
      // polarization J/psi
      normalize(_h_cTheta_Jpsi);
      pair<double,pair<double,double> > alpha = calcAlpha(_h_cTheta_Jpsi);
      Estimate1DPtr _h_alpha;
      book(_h_alpha,1,1,1);
      _h_alpha->bin(1).set(alpha.first, alpha.second);
      // polarization psi(2S)
      normalize(_h_cTheta_Psi_prime);
      alpha = calcAlpha(_h_cTheta_Psi_prime);
      book(_h_alpha,1,1,2);
      _h_alpha->bin(1).set(alpha.first, alpha.second);

      // J/psi as function of momentum
      book(_h_alpha,2,1,1);
      for (auto& b : _h2_cTheta_Jpsi->bins()) {
        normalize(b);
        alpha = calcAlpha(b);
        _h_alpha->bin(b.index()).set(alpha.first, alpha.second);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    // count of weights
    CounterPtr _weightSum;
    // histograms
    Histo1DPtr _h_Jpsi,_h_Psi_prime;
    Histo1DPtr _h_cTheta_Jpsi,_h_cTheta_Psi_prime;
    Histo1DGroupPtr _h2_cTheta_Jpsi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2002_I606309);

}
