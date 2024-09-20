// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_b0 -> Lambda_c+ pi-,K-
  class LHCB_2024_I2824757 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2024_I2824757);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==5122), "UFS" );
      for(unsigned int ip=0;ip<2;++ip) {
        for(unsigned int ix=0;ix<2;++ix) {
          for(unsigned int iy=0;iy<3;++iy) {
            string base = toString(ip)+"_"+toString(ix)+"_"+toString(iy);
            book(_h_cos1 [ip][ix][iy],"h_cos_1_" +base,20,-1,1);
            if(iy==2) continue;
            book(_h_cos2 [ip][ix][iy],"h_cos_2_" +base,20,-1,1);
            book(_p_cos12[ip][ix][iy],"p_cos_12_"+base,1,-1,1);
            book(_h_phi2 [ip][ix][iy],"h_phi_2_" +base,20,-M_PI,M_PI);
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_b0 baryons
      for (const Particle& Lamb : apply<UnstableParticles>(event, "UFS").particles()) {
        int sign = Lamb.pid()/5122;
        if(Lamb.children().size()!=2) continue;
        Particle baryon1,meson1;
        if ( Lamb.children()[0].pid()==sign*4122 &&
            (Lamb.children()[1].pid()==-sign*211 ||
             Lamb.children()[1].pid()==-sign*321)) {
          baryon1 = Lamb.children()[0];
          meson1  = Lamb.children()[1];
        }
        else if ( Lamb.children()[1].pid()==sign*4122 &&
                 (Lamb.children()[0].pid()==-sign*211 ||
                  Lamb.children()[0].pid()==-sign*321)) {
          baryon1 = Lamb.children()[1];
          meson1  = Lamb.children()[0];
        }
        else  continue;
        if(baryon1.children().size()!=2) continue;
        unsigned int ib = meson1.abspid()==211 ? 0 : 1;
        Particle baryon2,meson2;
        unsigned int ic=0;
        if (baryon1.children()[0].pid()== sign*3122 && baryon1.children()[1].pid()== sign*211) {
          baryon2 = baryon1.children()[0];
          meson2  = baryon1.children()[1];
                ic=0;
        }
        else if (baryon1.children()[1].pid()== sign*3122 && baryon1.children()[0].pid()== sign*211) {
          baryon2 = baryon1.children()[1];
          meson2  = baryon1.children()[0];
                ic=0;
        }
        else if (baryon1.children()[0].pid()== sign*3122 && baryon1.children()[1].pid()== sign*321) {
          baryon2 = baryon1.children()[0];
          meson2  = baryon1.children()[1];
                ic=1;
        }
        else if (baryon1.children()[1].pid()== sign*3122 && baryon1.children()[0].pid()== sign*321) {
          baryon2 = baryon1.children()[1];
          meson2  = baryon1.children()[0];
          ic=1;
        }
        else if( baryon1.children()[0].pid()== sign*2212 &&
                (baryon1.children()[1].pid()==-sign*311 ||
                 baryon1.children()[1].pid()== 310 )) {
          baryon2 = baryon1.children()[0];
          meson2  = baryon1.children()[1];
          ic=2;
        }
        else if( baryon1.children()[1].pid()== sign*2212 &&
                (baryon1.children()[1].pid()==-sign*311 ||
                 baryon1.children()[1].pid()== 310 )) {
          baryon2 = baryon1.children()[1];
          meson2  = baryon1.children()[0];
          ic=2;
        }
        else   continue;
        // deal with Kbar0 in Lambda_c+ -> Kbar0 p
        if (ic==2 && meson2.pid()!=310) {
          if (meson2.children().size()==1) {
            meson2=meson2.children()[0];
          }
          if (meson2.pid()!=310) continue;
        }
        // first boost to the Lambdab0 rest frame
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lamb.mom().betaVec());
        FourMomentum pbaryon1 = boost1.transform(baryon1.mom());
        FourMomentum pbaryon2 = boost1.transform(baryon2.mom());
        FourMomentum pmeson2  = boost1.transform(meson2.mom());
        // boost to lambda_c rest frame
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
        Vector3 axis1 = pbaryon1.p3().unit();
        FourMomentum pp = boost2.transform(pbaryon2);
        double cTheta1 = pp.p3().unit().dot(axis1);
        _h_cos1[(1-sign)/2][ib][ic]->fill(cTheta1);
        // that's it for Lambda_c+ -> p+ K_S0
        if (ic==2) continue;
        if (baryon2.children().size()!=2) continue;
        Particle baryon3,meson3;
        if (baryon2.children()[0].pid()== sign*2212 &&
            baryon2.children()[1].pid()==-sign*211) {
          baryon3 = baryon2.children()[0];
          meson3  = baryon2.children()[1];
        }
        else if (baryon2.children()[1].pid()== sign*2212 &&
          baryon2.children()[0].pid()==-sign*211) {
          baryon3 = baryon2.children()[1];
          meson3  = baryon2.children()[0];
        }
        else  continue;
        Vector3 axis2=pp.p3().unit();
        LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pp.betaVec());
        FourMomentum pbaryon3 = boost3.transform(boost2.transform(boost1.transform(baryon3.mom())));
        double cTheta2 = pbaryon3.p3().unit().dot(axis2);
        _h_cos2 [(1-sign)/2][ib][ic]->fill(cTheta2);
        _p_cos12[(1-sign)/2][ib][ic]->fill(cTheta1*cTheta2,cTheta1*cTheta2);
        Vector3 trans1 = axis1 - axis1.dot(axis2)*axis2;
        Vector3 trans2 = pbaryon3.p3()-pbaryon3.p3().dot(axis2)*pbaryon3.p3();
        double phi2 = atan2(trans1.cross(trans2).dot(axis2), trans1.dot(trans2));
        _h_phi2[(1-sign)/2][ib][ic]->fill(phi2);
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) const {
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

    vector<pair<double,double>> calcBetaGamma(Histo1DPtr hist) const {
      if(hist->numEntries()==0.) return vector<pair<double,double>>(2,make_pair(0.,0.));;
      double sum[7]={0.,0.,0.,0.,0.,0.,0.,};
      for (const auto& bin : hist->bins()) {
        double Oi = bin.sumW();
        if(Oi==0.) continue;
        double ai = (bin.xMax()-bin.xMin())/2./M_PI;
        double bi = M_PI/32.*(cos(bin.xMin())-cos(bin.xMax()));
        double ci =-M_PI/32.*(sin(bin.xMin())-sin(bin.xMax()));
        double Ei = bin.errW();
        sum[0] +=   bi*Oi/sqr(Ei);
        sum[1] +=   ci*Oi/sqr(Ei);
        sum[2] +=   ai*bi/sqr(Ei);
        sum[3] +=   ai*ci/sqr(Ei);
        sum[4] += sqr(bi)/sqr(Ei);
        sum[5] += sqr(ci)/sqr(Ei);
        sum[6] +=   bi*ci/sqr(Ei);
      }
      double gamma = (sum[0]/sum[4]-sum[2]/sum[4]-sum[1]/sum[6]+sum[3]/sum[6])/(sum[6]/sum[4]-sum[5]/sum[6]);
      double beta  = (sum[0]/sum[6]-sum[2]/sum[6]-sum[1]/sum[5]+sum[3]/sum[5])/(sum[4]/sum[6]-sum[6]/sum[5]);
      return {make_pair(beta ,1./sqrt(sum[4])),make_pair(gamma,1./sqrt(sum[5]))};
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_cos1);
      normalize(_h_cos2);
      normalize(_h_phi2);
      // compute values for specific modes
      pair<double,double> a_b[2][2][2], a_c[2][2][2], a_s[2][2][2],beta[2][2][2],gamma[2][2][2];
      for (unsigned int ip=0;ip<2;++ip) {
        for (unsigned int ib=0;ib<2;++ib) {
          for (unsigned int ic=0;ic<2;++ic) {
            pair<double,double> p1 = calcAlpha(_h_cos1[ip][ib][ic]);
            pair<double,double> p2 = calcAlpha(_h_cos2[ip][ib][ic]);
            pair<double,double> p3 = make_pair(9.*_p_cos12[ip][ib][ic]->bin(1).mean(2),
                                               9.*_p_cos12[ip][ib][ic]->bin(1).stdErr(2));
            double ferr = sqrt(sqr(p1.second/p1.first)+sqr(p2.second/p2.first)+sqr(p3.second/p3.first));
            // lamda decay first as sign well known
            double l3 = p2.first*p3.first/p1.first;
            double e3 = l3*ferr;
            l3 = sqrt(max(l3,0.));
            e3 = 0.5*e3;
            // sign ambiguity so fix alpha for lambda0 -> p pi to known sign
            if(ip==1) l3*=-1.;
            a_s[ip][ib][ic] = make_pair(l3,e3);
            // lambda_b decay
            double l1 = p1.first*p3.first/p2.first;
            double e1 = l1*ferr;
            l1 = sqrt(max(l1,0.));
            e1 = 0.5*e1;
            if(p3.first*l3<0.) l1 *= -1.;
            a_b[ip][ib][ic] = make_pair(l1,e1);
            // lambda_c decay
            double l2 = p1.first*p2.first/p3.first;
            double e2 = l2*ferr;
            l2 = sqrt(max(l2,0.));
            e2 = 0.5*e2;
            if(p2.first*l3<0.) l2 *= -1.;
            a_c[ip][ib][ic] = make_pair(l2,e2);
            vector<pair<double,double>> temp = calcBetaGamma(_h_phi2[ip][ib][ic]);
            for(auto & val : temp) {
              double error = sqrt(sqr(val.second/val.first)+sqr(p3.second/p3.first));
              val.first /=p3.first;
              val.second = val.first*error;
            }
            beta [ip][ib][ic] = temp[0];
            gamma[ip][ib][ic] = temp[1];
          }
        }
      }
      // perform averages
      // Lambda^0 -> p pi
      pair<double,double> alam[2];
      for(unsigned int ip=0;ip<2;++ip) {
        double sum1(0.),sum2(0.);
        for(unsigned int ib=0;ib<2;++ib) {
          for(unsigned int ic=0;ic<2;++ic) {
            sum1 += a_s[ip][ib][ic].first/sqr(a_s[ip][ib][ic].second);
            sum2 += 1./sqr(a_s[ip][ib][ic].second);
          }
        }
        alam[ip] = make_pair(sum1/sum2,sqrt(1./sum2));
        Estimate0DPtr tmp;
        book(tmp,1,6,1+ip);
        tmp->set(alam[ip].first,alam[ip].second);
      }
      Estimate0DPtr tmp;
      book(tmp,1,6,3);
      tmp->set(0.5*(alam[0].first-alam[1].first),
               0.5*sqrt(sqr(alam[0].second)+sqr(alam[1].second)));
      book(tmp,1,6,4);
      tmp->set((alam[0].first+alam[1].first)/(alam[0].first-alam[1].first),
               4.*(sqr(alam[0].first*alam[1].second)+sqr(alam[1].first*alam[0].second))/pow(alam[0].first-alam[1].first,4));
      // Lambda_b decays
      pair<double,double> ab[2][2];
      for(unsigned int ib=0;ib<2;++ib) {
        for(unsigned int ip=0;ip<2;++ip) {
          double sum1(0.),sum2(0.);
          for(unsigned int ic=0;ic<2;++ic) {
            sum1 += a_b[ip][ib][ic].first/sqr(a_b[ip][ib][ic].second);
            sum2 += 1./sqr(a_b[ip][ib][ic].second);
          }
          ab[ip][ib] = make_pair(sum1/sum2,sqrt(1./sum2));
          Estimate0DPtr tmp;
          book(tmp,1,1+ib,1+ip);
          tmp->set(ab[ip][ib].first,ab[ip][ib].second);
        }
        Estimate0DPtr tmp;
        book(tmp,1,1+ib,3);
        tmp->set(0.5*(ab[0][ib].first-ab[1][ib].first),
                 0.5*sqrt(sqr(ab[0][ib].second)+sqr(ab[1][ib].second)));
        book(tmp,1,1+ib,4);
        tmp->set((ab[0][ib].first+ab[1][ib].first)/(ab[0][ib].first-ab[1][ib].first),
                 4.*(sqr(ab[0][ib].first*ab[1][ib].second)+sqr(ab[1][ib].first*ab[0][ib].second))/pow(ab[0][ib].first-ab[1][ib].first,4));
      }
      // Lambda_c decays
      for(unsigned int ic=0;ic<2;++ic) {
        pair<double,double> ac[2],bc[2],gc[2];
        for(unsigned int ip=0;ip<2;++ip) {
          double sum1(0.),sum2(0.);
          for(unsigned int ib=0;ib<2;++ib) {
            sum1 += a_c[ip][ib][ic].first/sqr(a_c[ip][ib][ic].second);
            sum2 += 1./sqr(a_c[ip][ib][ic].second);
          }
          ac[ip] = make_pair(sum1/sum2,sqrt(1./sum2));
          Estimate0DPtr tmp;
          book(tmp,1,3+ic,1+ip);
          tmp->set(ac[ip].first,ac[ip].second);
          sum1=sum2=0.;
          for(unsigned int ib=0;ib<2;++ib) {
            sum1 += beta[ip][ib][ic].first/sqr(beta[ip][ib][ic].second);
            sum2 += 1./sqr(beta[ip][ib][ic].second);
          }
          bc[ip] = make_pair(sum1/sum2,sqrt(1./sum2));
          book(tmp,2,1+ip,1+ic);
          tmp->set(bc[ip].first,bc[ip].second);
          sum1=sum2=0.;
          for(unsigned int ib=0;ib<2;++ib) {
            sum1 += gamma[ip][ib][ic].first/sqr(gamma[ip][ib][ic].second);
            sum2 += 1./sqr(gamma[ip][ib][ic].second);
          }
          gc[ip] = make_pair(sum1/sum2,sqrt(1./sum2));
          book(tmp,2,3+ip,1+ic);
          tmp->set(gc[ip].first,gc[ip].second);





        }
        Estimate0DPtr tmp;
        book(tmp,1,3+ic,3);
        tmp->set(0.5*(ac[0].first-ac[1].first),
                 0.5*sqrt(sqr(ac[0].second)+sqr(ac[1].second)));
        book(tmp,1,3+ic,4);
        tmp->set((ac[0].first+ac[1].first)/(ac[0].first-ac[1].first),
                 4.*(sqr(ac[0].first*ac[1].second)+sqr(ac[1].first*ac[0].second))/pow(ac[0].first-ac[1].first,4));
      }
      // finally lambda_c -> p KS0
      pair<double,double> ap[2];
      for(unsigned int ip=0;ip<2;++ip) {
        double sum1(0.),sum2(0.);
        for(unsigned int ib=0;ib<2;++ib) {
          pair<double,double> p1 = calcAlpha(_h_cos1[ip][ib][2]);
          double val = p1.first/ab[ip][ib].first;
          double err = sqr(val)*(sqr(p1.second/p1.first) + sqr(ab[ip][ib].second/ab[ip][ib].first));
          sum1 += val/err;
          sum2 += 1./err;
        }
        ap[ip] = make_pair(sum1/sum2,sqrt(1./sum2));
        Estimate0DPtr tmp;
        book(tmp,1,5,1+ip);
        tmp->set(ap[ip].first,ap[ip].second);
      }
      book(tmp,1,5,3);
      tmp->set(0.5*(ap[0].first-ap[1].first),
               0.5*sqrt(sqr(ap[0].second)+sqr(ap[1].second)));
      book(tmp,1,5,4);
      tmp->set((ap[0].first+ap[1].first)/(ap[0].first-ap[1].first),
               4.*(sqr(ap[0].first*ap[1].second)+sqr(ap[1].first*ap[0].second))/pow(ap[0].first-ap[1].first,4));
    }
    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cos1[2][2][3],_h_cos2[2][2][2],_h_phi2[2][2][2];
    Profile1DPtr _p_cos12[2][2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2024_I2824757);

}
