// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > Lambda0 Lambdabar0
  class BESIII_2023_I2637702 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2637702);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::abspid==3122), "UFS");
      declare(FinalState(), "FS");

      // Book histograms
      book(_h_T2, "TMP/T2",20,-1.,1.);
      book(_h_T3, "TMP/T3",20,-1.,1.);
      book(_h_cThetaL,"cThetaL",20,-1.,1.);
      book(_wsum,"TMP/wsum");

      vector<string> energies({"3.68", "3.683", "3.684", "3.685", "3.687", "3.691", "3.71"});
      for(const string& en : energies) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(!inRange(sqrtS(),3.68,3.71) && _ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren (const Particle& p,map<long,int> & nRes, int &ncount) {
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
      // loop over lambda0 baryons
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      Particle Lambda,LamBar;
      bool matched(false);
      for (const Particle& p : ufs.particles(Cuts::abspid==3122)) {
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
                Lambda = p;
                LamBar = p2;
              }
              else {
                Lambda = p2;
                LamBar = p;
              }
       	      break;
       	    }
       	  }
       	}
      	if (matched) break;
      }
      if (!matched) vetoEvent;
      Particle proton;
      matched = false;
      for (const Particle& p : Lambda.children()) {
        if (p.pid()==2212) {
          matched=true;
          proton=p;
        }
        else if (p.pid()==PID::PHOTON) {
          vetoEvent;
        }
      }
      if (!matched) vetoEvent;
      Particle baryon;
      matched = false;
      for (const Particle & p : LamBar.children()) {
        if (p.pid()==-2212) {
          baryon=p;
          matched=true;
        }
        else if (p.pid()==PID::PHOTON) {
          vetoEvent;
        }
      }
      if (!matched) vetoEvent;
      // boost to the Lambda rest frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambda.mom().betaVec());
      Vector3 e1z = Lambda.mom().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      Vector3 axis1 = boost1.transform(proton.mom()).p3().unit();
      double n1x(e1x.dot(axis1)),n1y(e1y.dot(axis1)),n1z(e1z.dot(axis1));
      // boost to the Lambda bar
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(LamBar.mom().betaVec());
      Vector3 axis2 = boost2.transform(baryon.mom()).p3().unit();
      double n2x(e1x.dot(axis2)),n2z(e1z.dot(axis2));
      double cosL = axis.dot(Lambda.mom().p3().unit());
      double sinL = sqrt(1.-sqr(cosL));
      double T2 = -sinL*cosL*(n1x*n2z+n1z*n2x);
      double T3 = -sinL*cosL*n1y;
      _h_T2->fill(cosL,T2);
      _h_T3->fill(cosL,T3);
      _wsum->fill();
      _h_cThetaL->fill(cosL);
    }

    pair<double,pair<double,double> > calcAlpha0(Histo1DPtr hist) {
      if (hist->numEntries()==0.)  return make_pair(0.,make_pair(0.,0.));
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

    pair<double,double> calcCoeff(unsigned int imode, Histo1DPtr hist) {
      if (hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins()) {
        double Oi = bin.sumW();
        if (Oi==0.) continue;
        double ai(0.),bi(0.);
        if (imode==0) {
          bi = (pow(1.-sqr(bin.xMin()),1.5) - pow(1.-sqr(bin.xMax()),1.5))/3.;
        }
        else if(imode>=2 && imode<=4) {
          bi = (pow(bin.xMin(),3)*( -5. + 3.*sqr(bin.xMin())) + pow(bin.xMax(),3)*(5. - 3.*sqr(bin.xMax())))/15.;
        }
        else {
          assert(false);
        }
        double Ei = bin.errW();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double aLambda = 0.754;
      normalize(_h_cThetaL);
      scale(_h_T2,1./ *_wsum);
      scale(_h_T3,1./ *_wsum);
      pair<double,pair<double,double> > alpha0 = calcAlpha0(_h_cThetaL);
      pair<double,pair<double,double> > R;
      double tau = sqr(sqrtS()/(2*1.115683));
      R.first = sqrt(tau*(1-alpha0.first)/(1+alpha0.first));
      R.second.first = R.second.second = R.first/(1.-sqr(alpha0.first));

      pair<double,double> c_T2 = calcCoeff(2,_h_T2);
      pair<double,double> c_T3 = calcCoeff(3,_h_T3);

      double sDelta = (-2.*(3. + alpha0.first)*c_T3.first)/(-aLambda*sqrt(1 - sqr(alpha0.first)));
      double cDelta = (-3*(3 + alpha0.first)*c_T2.first)/(-sqr(aLambda)*sqrt(1 - sqr(alpha0.first)));

      pair<double,pair<double,double> > Delta;
      Delta.first = asin(sDelta);
      if(cDelta<0.) Delta.first = M_PI-Delta.first;
      Delta.second.first = (-4*(sqr(c_T3.second)*sqr(1. + alpha0.first)*
				sqr(1 + alpha0.first)*sqr(3.+alpha0.first) +
				sqr(alpha0.second.first)*sqr(c_T3.first)*sqr(1 + 3*alpha0.first)))
	/(sqr(1.-sqr(alpha0.first))*(4*sqr(c_T3.first)*sqr(3 + alpha0.first) + sqr(aLambda)*
				     (-1 + sqr(alpha0.first))));
      Delta.second.second = (-4*(sqr(c_T3.second)*sqr(1. + alpha0.first)*
				 sqr(1 + alpha0.first)*sqr(3.+alpha0.first) +
				 sqr(alpha0.second.second)*sqr(c_T3.first)*sqr(1 + 3*alpha0.first)))
	/(sqr(1.-sqr(alpha0.first))*(4*sqr(c_T3.first)*sqr(3 + alpha0.first) + sqr(aLambda)*
				     (-1 + sqr(alpha0.first))));
      Delta.first         *=180./M_PI;
      Delta.second.first  *=180./M_PI;
      Delta.second.second *=180./M_PI;
      for(unsigned int iy=0;iy<3;++iy) {
	double val;
	pair<double,double> err;
	if(iy==0) {
	  val = alpha0.first;
	  err = alpha0.second;
	}
	else if(iy==1) {
	  val = Delta.first ;
	  err = Delta.second;
	}
	else {
	  val = R.first;
	  err = R.second;
	}
        Estimate1DPtr tmp1;
        book(tmp1,1,1,1+iy);
        tmp1->bin(1).set(val,err);
        BinnedEstimatePtr<string> tmp2;
        book(tmp2,2,1,1+iy);
        tmp2->binAt(_ecms).set(val,err);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_T2,_h_T3;
    Histo1DPtr _h_cThetaL;
    CounterPtr _wsum;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2637702);

}
