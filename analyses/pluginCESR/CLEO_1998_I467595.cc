// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D* helicity in e+e- at 10.5 GeV
  class CLEO_1998_I467595 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1998_I467595);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_ctheta, {0.25, 0.45, 0.55, 0.65, 0.75, 0.85, 1.});
      for (auto& b : _h_ctheta->bins()) {
        book(b, 5, 1, b.index());
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.01026));

      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==413)) {
        if (p.children().size()!=2) continue;
        int sign = p.pid()/413;
        Particle D0;
        if(p.children()[0].pid()==sign*421 && p.children()[1].pid()==sign*211) {
          D0 = p.children()[0];
        }
        else if(p.children()[1].pid()==sign*421 && p.children()[0].pid()==sign*211) {
          D0 = p.children()[1];
        }
        else {
          continue;
        }
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        double xp = (p.momentum().p3().mod()+p.momentum().t())/(Pmax+Emax);
        Vector3 e1z = p.momentum().p3().unit();
        Vector3 axis1 = boost.transform(D0.momentum()).p3().unit();
        double ctheta = e1z.dot(axis1);
        _h_ctheta->fill(xp,ctheta);
      }
    }

    pair<double,double> calcRho(const Histo1DPtr& hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins() ) {
        double Oi = bin.sumW();
        if(Oi==0.) continue;
        double ai = 0.25*(bin.xMax()*(3.-sqr(bin.xMax())) - bin.xMin()*(3.-sqr(bin.xMin())));
        double bi = 0.75*(bin.xMin()*(1.-sqr(bin.xMin())) - bin.xMax()*(1.-sqr(bin.xMax())));
        double Ei = bin.errW();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
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
      if(dis>0.) {
        dis = sqrt(dis);
        return make_pair(alpha,make_pair(0.5*(-bb+dis)/aa,-0.5*(-bb-dis)/aa));
      }
      else {
        return make_pair(alpha,make_pair(0.,0.));
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      vector<double> x = {0.25,0.45,0.55,0.65,0.75,0.85,1.};
      Estimate1DPtr h_alpha;
      book(h_alpha, 3,1,1);
      Estimate1DPtr h_rho;
      book(h_rho, 4,1,1);
      for (auto& b : _h_ctheta->bins()) {
        normalize(b);
        pair<double,double> rho00 = calcRho(b);
        h_rho->bin(b.index()).set(rho00.first, rho00.second);
        pair<double,pair<double,double> > alpha = calcAlpha(b);
        h_alpha->bin(b.index()).set(alpha.first, alpha.second);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_ctheta;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_1998_I467595);


}
