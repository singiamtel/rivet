// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > Lambda_c+ Lambda_c-
  class BESIII_2023_I2677290 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2677290);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::abspid==4122), "UFS");
      declare(FinalState(), "FS");
      // histograms
      book(_h_cThetaL,"cThetaL",20,-1.,1.);
      book(_wsum,"TMP/wsum");

      vector<string> energies({"4.6119", "4.628", "4.6409", "4.6612", "4.6819", "4.6988",
          "4.7397", "4.75", "4.7805", "4.8431", "4.918", "4.9509"});
      for(const string& en : energies) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int & ncount) {
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
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particle Lambda, LamBar;
      bool matched(false);
      for (const Particle& p :  ufs.particles()) {
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
      const double cosL = axis.dot(Lambda.mom().p3().unit());
      _wsum->fill();
      _h_cThetaL->fill(cosL);
    }

    pair<double,pair<double,double> > calcAlpha0(Histo1DPtr hist) {
      if (hist->numEntries()==0.)  return make_pair(0.,make_pair(0.,0.));
      const double d = 3./(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      const double c = 3.*(hist->xMax()-hist->xMin())/(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double sum1(0.),sum2(0.),sum3(0.),sum4(0.),sum5(0.);
      for (const auto& bin : hist->bins()) {
       	const double Oi = bin.sumW();
        if (Oi==0.) continue;
        const double a =  d*(bin.xMax() - bin.xMin());
        const double b = d/3.*(pow(bin.xMax(),3) - pow(bin.xMin(),3));
       	const double Ei = bin.errW();
        sum1 +=   a*Oi/sqr(Ei);
        sum2 +=   b*Oi/sqr(Ei);
        sum3 += sqr(a)/sqr(Ei);
        sum4 += sqr(b)/sqr(Ei);
        sum5 +=    a*b/sqr(Ei);
      }
      // calculate alpha
      const double alpha = (-c*sum1 + sqr(c)*sum2 + sum3 - c*sum5)/(sum1 - c*sum2 + c*sum4 - sum5);
      // and error
      const double cc = -pow((sum3 + sqr(c)*sum4 - 2*c*sum5),3);
      const double bb = -2*sqr(sum3 + sqr(c)*sum4 - 2*c*sum5)*(sum1 - c*sum2 + c*sum4 - sum5);
      const double aa =  sqr(sum1 - c*sum2 + c*sum4 - sum5)*(-sum3 - sqr(c)*sum4 + sqr(sum1 - c*sum2 + c*sum4 - sum5) + 2*c*sum5);
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
      // storage of the values to fill histos
      const double mn    = 2.28646;
      const double tau   = 4.*sqr(mn/sqrtS());
      const double beta  = sqrt(1.-tau);
      const double alpha = 7.2973525693e-3;
      const double GeV2pb = 0.3893793721e9;
      scale(_wsum, crossSection()/ sumOfWeights()/picobarn);
      // prefactor and sigma
      double sigma0 = 4.*M_PI*sqr(alpha/sqrtS())*beta*GeV2pb;
      // calculate alpha0
      pair<double,pair<double,double> > alpha0 = calcAlpha0(_h_cThetaL);
      // Geff
      pair<double,double> Geff = make_pair(1e2*sqrt(3.*_wsum->val()/(sigma0*(1 + 0.5*tau))),
                                           1e2*1.5*_wsum->val()/(sigma0*(1 + 0.5*tau)));
      // GM
      double GMv = 1e2*sqrt(6.*((1+alpha0.first)*_wsum->val())/((3+alpha0.first)*sigma0));
      double GMe1 = 1e2*sqrt((3*(sqr(_wsum->err())*sqr(1 + alpha0.first)*sqr(3 + alpha0.first) + 4*sqr(alpha0.second.first)*sqr(_wsum->val())))/
                             (2.*(1 + alpha0.first)*pow(3 + alpha0.first,3)*_wsum->val()*sigma0));
      double GMe2 = 1e2*sqrt((3*(sqr(_wsum->err())*sqr(1 + alpha0.first)*sqr(3 + alpha0.first) + 4*sqr(alpha0.second.second)*sqr(_wsum->val())))/
                             (2.*(1 + alpha0.first)*pow(3 + alpha0.first,3)*_wsum->val()*sigma0));
      pair<double,pair<double,double> > GM = make_pair(std::move(GMv), make_pair(std::move(GMe1), std::move(GMe2)));
      // ratio
      pair<double,pair<double,double>> R;
      R.first = sqrt((1 - alpha0.first)/(tau + alpha0.first*tau));
      R.second.first  = R.first*alpha0.second.first /(1.-sqr(alpha0.first));
      R.second.second = R.first*alpha0.second.second/(1.-sqr(alpha0.first));
      for(unsigned int ix=1; ix<6; ++ix) {
        double val;
        pair<double,double> err;
        if (ix==1) {
          val = _wsum->val();
          err = make_pair(_wsum->err(),_wsum->err());
        }
        else if (ix==2) {
          val = Geff.first;
          err = make_pair(Geff.second,Geff.second);
        }
        else if (ix==3) {
          val = alpha0.first;
          err = alpha0.second;
        }
        else if (ix==4) {
          val = R.first;
          err = R.second;
        }
        else {
          val = GM.first;
          err = GM.second;
        }
        BinnedEstimatePtr<string> tmp;
        book(tmp, 1, 1, ix);
        tmp->binAt(_ecms).set(val,err);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cThetaL;
    CounterPtr _wsum;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2677290);

}
