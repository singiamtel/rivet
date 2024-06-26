// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> Xi0 Xibar0
  class BESIII_2023_I2634735 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2634735);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");
      // Book histograms
      book(_h_T1, "TMP/T1",20,-1.,1.);
      book(_h_T2, "TMP/T2",20,-1.,1.);
      book(_h_T3, "TMP/T3",20,-1.,1.);
      book(_h_T4, "TMP/T4",20,-1.,1.);
      book(_h_T5, "TMP/T5",20,-1.,1.);
      book(_h_cTheta,"TMP/cTheta",20,-1.,1.);
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_cProton[ix],"TMP/cProton_"+toString(ix+1), 20, -1., 1.);
      }
      book(_wsum,"TMP/wsum");
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int& ncount) {
      for (const Particle& child : p.children()) {
        if (child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.pid()>0) {
      	axis = beams.first .mom().p3().unit();
      }
      else {
      	axis = beams.second.mom().p3().unit();
      }
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // loop over Xi baryons
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particle Xi,XiBar;
      bool matched(false);
      for (const Particle& p : ufs.particles(Cuts::abspid==3322)) {
        if (p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        matched=false;
        // check for antiparticle
        for (const Particle& p2 :  ufs.particles(Cuts::pid==-p.pid())) {
          if (p2.children().empty()) continue;
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2==0) {
            matched = true;
            for (const auto& val : nRes2) {
              if (val.second!=0) {
                matched = false;
                break;
              }
            }
            // found baryon and antibaryon
            if (matched) {
              if (p.pid()>0) {
                Xi    = p;
                XiBar = p2;
              }
              else {
                Xi    = p2;
                XiBar = p;
              }
              break;
            }
          }
        }
        if (matched) break;
      }
      if (!matched) vetoEvent;
      // find the lambda and antilambda
      Particle Lambda,LamBar;
      if (Xi.children()[0].pid() ==3122 && Xi.children()[1].pid()==111) {
      	Lambda = Xi.children()[0];
      }
      else if (Xi.children()[1].pid() ==3122  && Xi.children()[0].pid()==111) {
      	Lambda = Xi.children()[1];
      }
      else {
        vetoEvent;
      }
      if (XiBar.children()[0].pid() ==-3122  && XiBar.children()[1].pid()==111) {
      	LamBar = XiBar.children()[0];
      }
      else if (XiBar.children()[1].pid() ==-3122 && XiBar.children()[0].pid()==111) {
      	LamBar = XiBar.children()[1];
      }
      else {
        vetoEvent;
      }
      // boost to the Xi rest frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Xi.mom().betaVec());
      Vector3 e1z = Xi.mom().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      FourMomentum pLambda = boost1.transform(Lambda.mom());
      Vector3 axis1 = pLambda.p3().unit();
      double n1x(e1x.dot(axis1)),n1y(e1y.dot(axis1)),n1z(e1z.dot(axis1));
      Particle proton;
      if (Lambda.children().size()!=2) vetoEvent;
      if (Lambda.children()[0].pid()== 2212 && Lambda.children()[1].pid()==-211) {
        proton = Lambda.children()[0];
      }
      else if (Lambda.children()[1].pid()== 2212 && Lambda.children()[0].pid()==-211) {
        proton = Lambda.children()[1];
      }
      else {
        vetoEvent;
      }
      LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pLambda.betaVec());
      FourMomentum pProton = boost3.transform(boost1.transform(proton.mom()));
      const double cProton = pProton.p3().unit().dot(axis1);
      _h_cProton[0]->fill(cProton);
      // boost to the Xi bar rest frame
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(XiBar.mom().betaVec());
      FourMomentum pLamBar = boost2.transform(LamBar.mom());
      Vector3 axis2 = pLamBar.p3().unit();
      double n2x(e1x.dot(axis2)),n2y(e1y.dot(axis2)),n2z(e1z.dot(axis2));
      double cosX = axis.dot(Xi.mom().p3().unit());
      double sinX = sqrt(1.-sqr(cosX));
      Particle pbar;
      if (LamBar.children().size()!=2) vetoEvent;
      if (LamBar.children()[0].pid()==-2212 && LamBar.children()[1].pid()== 211) {
        pbar = LamBar.children()[0];
      }
      else if (LamBar.children()[1].pid()==-2212 && LamBar.children()[0].pid()== 211) {
        pbar = LamBar.children()[1];
      }
      else {
        vetoEvent;
      }
      LorentzTransform boost4 = LorentzTransform::mkFrameTransformFromBeta(pLamBar.betaVec());
      FourMomentum pPbar = boost4.transform(boost2.transform(pbar.mom()));
      double cPbar = pPbar.p3().unit().dot(axis2);
      _h_cProton[1]->fill(cPbar);
      // moments
      double T1 = sqr(sinX)*n1x*n2x+sqr(cosX)*n1z*n2z;
      double T2 = -sinX*cosX*(n1x*n2z+n1z*n2x);
      double T3 = -sinX*cosX*n1y;
      double T4 = -sinX*cosX*n2y;
      double T5 = n1z*n2z-sqr(sinX)*n1y*n2y;
      _h_T1->fill(cosX,T1);
      _h_T2->fill(cosX,T2);
      _h_T3->fill(cosX,T3);
      _h_T4->fill(cosX,T4);
      _h_T5->fill(cosX,T5);
      _h_cTheta->fill(cosX);
      _wsum->fill();
    }

    pair<double,pair<double,double> > calcAlpha0(Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,make_pair(0.,0.));
      double d = 3./(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double c = 3.*(hist->xMax()-hist->xMin())/(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double sum1(0.),sum2(0.),sum3(0.),sum4(0.),sum5(0.);
      for (const auto& bin : hist->bins() ) {
       	double Oi = bin.sumW();
        if (Oi==0.) continue;
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

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.), sum2(0.);
      for (const auto& bin : hist->bins() ) {
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

    pair<double,double> calcCoeff(unsigned int imode, Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins()) {
        double Oi = bin.sumW();
        if (Oi==0.) continue;
        double ai(0.), bi(0.);
        if (imode==0) {
          bi = (pow(1.-sqr(bin.xMin()),1.5) - pow(1.-sqr(bin.xMax()),1.5))/3.0;
        }
        else if (imode>=2 && imode<=4) {
          bi = (pow(bin.xMin(),3)*(-5. + 3.*sqr(bin.xMin())) + pow(bin.xMax(),3)*(5. - 3.*sqr(bin.xMax())))/15.;
        }
        else {
          assert(false);
        }
        const double Ei = bin.errW();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double aLambda = 0.754;
      normalize(_h_cTheta);
      normalize(_h_cProton);
      scale(_h_T1, 1.0/ *_wsum);
      scale(_h_T2, 1.0/ *_wsum);
      scale(_h_T3, 1.0/ *_wsum);
      scale(_h_T4, 1.0/ *_wsum);
      scale(_h_T5, 1.0/ *_wsum);
      // calculate alpha0
      pair<double,pair<double,double> > alpha0 = calcAlpha0(_h_cTheta);
      Estimate0DPtr _h_alpha0;
      book(_h_alpha0,1,1,1);
      _h_alpha0->set(alpha0.first,  alpha0.second);
      double s2 = -1. + sqr(alpha0.first);
      double s3 = 3 + alpha0.first;
      double s1 = sqr(s3);
      // alpha parameters
      pair<double,double> alpha[2];
      for(unsigned int ix=0;ix<2;++ix) {
	Estimate0DPtr _h_alpha;
	book(_h_alpha,1,1,3+2*ix);
	alpha[ix] = calcAlpha(_h_cProton[ix]);
	alpha[ix].first  /= aLambda;
	if(ix==1) alpha[ix].first *=-1;
	alpha[ix].second /= aLambda;
	_h_alpha->set(alpha[ix].first, alpha[ix].second);
      }
      // now for Delta
      pair<double,double> c_T2 = calcCoeff(2,_h_T2);
      pair<double,double> c_T3 = calcCoeff(3,_h_T3);
      pair<double,double> c_T4 = calcCoeff(4,_h_T4);
      double s4 = sqr(c_T2.first);
      double s5 = sqr(c_T3.first);
      double s6 = sqr(c_T4.first);
      double disc = s1*s5*s6*(-9.*s2*s4 + 4.*s1*s5*s6);
      if (disc<0.)  return;
      disc = sqrt(disc);
      double sDelta = (-2.*(3. + alpha0.first)*c_T3.first)/(alpha[0].first*sqrt(1 - sqr(alpha0.first)));
      double cDelta = (-3*(3 + alpha0.first)*c_T2.first)/(alpha[0].first*alpha[1].first*sqrt(1 - sqr(alpha0.first)));
      double Delta = asin(sDelta);
      if(cDelta<0.) Delta = M_PI-Delta;
      double ds_P = (-9*c_T2.first*((-1 + alpha0.first)*(1 + alpha0.first) *
                    (3 + alpha0.first)*c_T3.first*c_T4.first*c_T2.second +
                    c_T2.first*c_T4.first*(c_T3.first*(alpha0.second.first +
                    3*alpha0.first*alpha0.second.first) -(-1 + alpha0.first) *
                    (1 + alpha0.first)*(3 + alpha0.first)*c_T3.second) -
                    (-1 + alpha0.first)*(1 + alpha0.first) *
                    (3 + alpha0.first)*c_T2.first*c_T3.first*c_T4.second)*disc) /
                    (pow(1 - pow(alpha0.first,2),1.5)*pow(c_T4.first,3) *
                    pow(-((disc + 2*s1*s5*s6) / (s2*s6)),1.5)*(-9*s2*s4 + 4*s1*s5*s6));
      double ds_M = (-9*c_T2.first*((-1 + alpha0.first)*(1 + alpha0.first) *
                    (3 + alpha0.first)*c_T3.first*c_T4.first*c_T2.second +
                    c_T2.first*c_T4.first*(c_T3.first*(alpha0.second.second +
                    3*alpha0.first*alpha0.second.second) - (-1 + alpha0.first) *
                    (1 + alpha0.first)*(3 + alpha0.first)*c_T3.second) -
                    (-1 + alpha0.first)*(1 + alpha0.first) *
                    (3 + alpha0.first)*c_T2.first*c_T3.first * c_T4.second)*disc) /
                    (pow(1 - pow(alpha0.first,2),1.5)*pow(c_T4.first,3) *
                    pow(-((disc + 2*s1*s5*s6) /(s2*s6)),1.5)*(-9*s2*s4 + 4*s1*s5*s6));
      ds_P /= sqrt(1.-sqr(sDelta));
      ds_M /= sqrt(1.-sqr(sDelta));
      Estimate0DPtr _h_sin;
      book(_h_sin, 1, 1, 2);
      _h_sin->set(Delta, make_pair( -ds_P, -ds_M));
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_T1,_h_T2,_h_T3,_h_T4,_h_T5;
    Histo1DPtr _h_cTheta,_h_cProton[2];
    CounterPtr _wsum;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2634735);

}
