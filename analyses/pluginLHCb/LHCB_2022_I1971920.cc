// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambdab0 -> Lambda0 gamma
  class LHCB_2022_I1971920 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2022_I1971920);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==5122), "UFS" );
      // histos
      book(_h_ctheta_gamma, "TMP/ctheta_gamma", 20, -1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_b baryons
      for (const Particle& lamB : apply<UnstableParticles>(event, "UFS").particles()) {
       	int sign = lamB.pid()/5122;
        // Lambda_b -> Lambda0 gamma
      	if (lamB.children().size()!=2) continue;
       	Particle lambda,gamma;
        if (lamB.children()[1].pid()==sign*3122 &&
           lamB.children()[0].pid()==22) {
          lambda = lamB.children()[1];
          gamma  = lamB.children()[0];
        }
        else if(lamB.children()[0].pid()==sign*3122 &&
          lamB.children()[1].pid()==22) {
          lambda = lamB.children()[0];
          gamma  = lamB.children()[1];
        }
        else continue;
        // Lambda0 -> p pi+
      	if (lambda.children().size()!=2) continue;
       	Particle proton, pion;
        if (lambda.children()[0].pid()== sign*2212 &&
           lambda.children()[1].pid()==-sign*211) {
          proton = lambda.children()[0];
          pion   = lambda.children()[1];
        }
        else if (lambda.children()[1].pid()== sign*2212 &&
          lambda.children()[0].pid()==-sign*211) {
          proton = lambda.children()[1];
          pion   = lambda.children()[0];
        }
        else  continue;
      	// first boost to the lamB rest frame
      	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(lamB.momentum().betaVec());
      	FourMomentum plambda = boost1.transform(lambda.mom());
      	FourMomentum pproton = boost1.transform(proton.mom());
      	// to lambda rest frame
      	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(plambda.betaVec());
      	Vector3 axis = plambda.p3().unit();
      	FourMomentum pp = boost2.transform(pproton);
      	// calculate angle
      	const double cTheta = pp.p3().unit().dot(axis);
        _h_ctheta_gamma->fill(cTheta);
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins()) {
        double Oi = bin.sumW();
        if (Oi==0.) continue;
        double ai = 0.5*(bin.xMax()-bin.xMin());
        double bi = 0.5*ai*(bin.xMax()+bin.xMin());
        double Ei = bin.errW();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_ctheta_gamma);
      pair<double,double> alpha = calcAlpha(_h_ctheta_gamma);
      const double aLam = 0.754;
      alpha.first  /= aLam;
      alpha.second /= aLam;
      Estimate0DPtr tmp;
      book(tmp,1,1,1);
      tmp->set(alpha.first,alpha.second);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_ctheta_gamma;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2022_I1971920);

}
