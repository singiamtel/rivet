// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_cJ -> phi phi
  class BESIII_2023_I2627838 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2627838);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(Cuts::abspid==100443), "UFS");
      // counters
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_n[ix],"TMP/n_"+toString(ix+1));
        for (unsigned int iy=0; iy<ix+1; ++iy)
          book(_m[ix][iy],"TMP/m_"+toString(ix+1)+"_"+toString(iy+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& psi : apply<UnstableParticles>(event, "UFS").particles()) {
        if (psi.children().size()!=2) continue;
        // find chi_c gamma decay and type of chi_c
        Particle chi;
        if (psi.children()[0].pid()==PID::GAMMA) {
          if (psi.children()[1].pid()==10441 ||
              psi.children()[1].pid()==20443 ||
              psi.children()[1].pid()==445 ) {
            chi = psi.children()[1];
          }
        }
        else if (psi.children()[1].pid()==PID::GAMMA) {
          if (psi.children()[0].pid()==10441 ||
              psi.children()[0].pid()==20443 ||
              psi.children()[0].pid()==445) {
	          chi = psi.children()[0];
          }
        }
        else {
          continue;
        }
        unsigned int iloc=0;
        if      (chi.pid()==10441) iloc=0;
        else if (chi.pid()==20443) iloc=1;
        else if (chi.pid()==  445) iloc=2;
        else {
          continue;
        }
        // require chi_c -> phi phi
        if (chi.children().size()!=2) continue;
        if (chi.children()[0].pid()!=333 ||
            chi.children()[0].pid()!=333) continue;
        bool found = true;
        Particle Km[2],Kp[2];
        for (unsigned int ix=0; ix<2; ++ix) {
          // required K+K- decay
          if (chi.children()[ix].children().size()!=2) {
            found = false;
            break;
          }
          if (chi.children()[ix].children()[0].pid()!=-chi.children()[ix].children()[1].pid() ||
              chi.children()[ix].children()[0].abspid()!=321) {
            found = false;
            break;
          }
          if (chi.children()[ix].children()[0].pid()>0) {
            Kp[ix] = chi.children()[ix].children()[0];
            Km[ix] = chi.children()[ix].children()[1];
          }
          else {
            Kp[ix] = chi.children()[ix].children()[1];
            Km[ix] = chi.children()[ix].children()[0];
          }
        }
        if (!found) continue;
        // fill count
        _n[iloc]->fill();
        // boost to psi(2S) frame
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(psi.mom().betaVec());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(chi.mom().betaVec());
        Vector3 axis1 = boost1.transform(chi.mom()).p3().unit();
        double cTheta[3]={0.,0.,0.};
        for (unsigned int ix=0; ix<2; ++ix) {
          FourMomentum pPhi = boost2.transform(boost1.transform(chi.children()[ix].mom()));
          Vector3 axis2 = pPhi.p3().unit();
          if (ix==0) cTheta[0] = axis1.dot(axis2);
          LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
          FourMomentum pK = boost3.transform(boost2.transform(boost1.transform(Kp[ix].mom())));
          cTheta[ix+1] = axis2.dot(pK.p3().unit());
        }
        if (iloc==0) {
          _m[0][0]->fill(0.25*(3. - 5.*sqr(cTheta[1])));
        }
        else if (iloc==1) {
          _m[1][0]->fill(0.625*(1. + sqr(cTheta[1]) - 4.*sqr(cTheta[0])));
          _m[1][1]->fill(     -(3.                  -10.*sqr(cTheta[0])));
        }
        else {
          _m[2][0]->fill(-0.125*(25.*sqr(cTheta[1])*sqr(cTheta[2]) - 10*sqr(cTheta[1]) - 10*sqr(cTheta[2]) + 3.));
          _m[2][1]->fill( 0.05 *(8 - 140.*sqr(cTheta[1]) + 325.*sqr(cTheta[0])*sqr(cTheta[2])
                         + 250.*sqr(cTheta[0])*sqr(cTheta[1])*sqr(cTheta[2])));
          _m[2][2]->fill( 0.75 *(1. + 5.*sqr(cTheta[1]) - 25.* sqr(cTheta[0])*sqr(cTheta[2])));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<3; ++ix) {
        if (_n[ix]->numEntries()==0)  continue;
        scale(_m[ix], 1.0/ *_n[ix]);
        if (ix==0) {
          double x = _m[0][0]->val()/(1.-2.*_m[0][0]->val());
          pair<double,double> dx = make_pair(x-(-_m[0][0]->err() + _m[0][0]->val())/(1 + 2*_m[0][0]->err() - 2*_m[0][0]->val()),
                                            (_m[0][0]->err() + _m[0][0]->val())/(1 - 2*_m[0][0]->err() - 2*_m[0][0]->val())-x);
          double rx = sqrt(abs(x));
          dx.first  *= 0.5/rx;
          dx.second *= 0.5/rx;
          if (x<0.) rx*=-1.;
          Estimate0DPtr h_x;
          book(h_x, 1, 1, 1);
          h_x->set(rx, dx);
        }
        else if (ix==1) {
          double u1 = -4.*_m[1][0]->val()/(-1. + 4.*_m[1][0]->val() + _m[1][1]->val());
          double u2 = -4.*_m[1][1]->val()/(-1. + 4.*_m[1][0]->val() + _m[1][1]->val());
          double O1 = _m[1][0]->val(), DO1 = _m[1][0]->err();
          double O2 = _m[1][1]->val(), DO2 = _m[1][1]->err();
          double root1 = sqrt(sqr(-1 + 4*O1 + O2)*(sqr(DO2)*sqr(1 - 4*O1) + 16*sqr(DO1)*sqr(O2)));
          pair<double,double> du1 = make_pair(-4*DO1*DO2/(4*DO1*DO2*(-1 + 4*O1 + O2) - root1),
                                               4*DO1*DO2/(4*DO1*DO2*(-1 + 4*O1 + O2) + root1));
          double ru1 = sqrt(abs(u1));
          du1.first  *=0.5/ru1;
          du1.second *=0.5/ru1;
          if (u1<0.) ru1*=-1.;
          double root2 = sqrt((sqr(DO2)*sqr(O1) + sqr(DO1)*sqr(-1 + O2))*sqr(-1 + 4*O1 + O2));
          pair<double,double> du2 = make_pair(-4*DO1*DO2/(DO1*DO2*(-1 + 4*O1 + O2) - root2),
                                               4*DO1*DO2/(DO1*DO2*(-1 + 4*O1 + O2) + root2));
          double ru2 = sqrt(abs(u2));
          du2.first  *=0.5/ru2;
          du2.second *=0.5/ru2;
          if (u2<0.) ru2*=-1.;
          Estimate0DPtr h_u1;
          book(h_u1, 1, 2, 1);
          h_u1->set(ru1, du1);
          Estimate0DPtr h_u2;
          book(h_u2, 1, 2, 2);
          h_u2->set(ru2, du2);
        }
        else {
          double O1 = _m[2][0]->val(), DO1 = _m[2][0]->err();
          double O2 = _m[2][1]->val(), DO2 = _m[2][1]->err();
          double O4 = _m[2][2]->val(), DO4 = _m[2][2]->err();
          double w1 = -O1/(-1 + 4*O1 + 2*O2 + 2*O4);
          double w2 = -O2/(-1 + 4*O1 + 2*O2 + 2*O4);
          double w4 = -O4/(-1 + 4*O1 + 2*O2 + 2*O4);
          double root1 = sqrt((16*sqr(DO1)*sqr(DO4)*sqr(O2) + sqr(DO2)*(sqr(DO4) *
                         sqr(1 - 4*O1) + 16*sqr(DO1)*sqr(O4)))/sqr(-1 + 4*O1 + 2*O2 + 2*O4));
          pair<double,double> dw1 = make_pair( DO1*DO2*DO4*(4*DO1*DO2*DO4 + (-1 + 4*O1 + 2*O2 + 2*O4)*root1) /
                                               ((-1 + 4*O1 + 2*O2 + 2*O4)*(16*sqr(DO1)*sqr(DO4)*sqr(O2) +
                                               sqr(DO2)*(sqr(DO4)*(-16*sqr(DO1) + sqr(1 - 4*O1)) + 16*sqr(DO1)*sqr(O4)))),
                                               DO1*DO2*DO4*(-4*DO1*DO2*DO4 + (-1 + 4*O1 + 2*O2 + 2*O4)*root1) /
                                               ((-1 + 4*O1 + 2*O2 + 2*O4)*(16*sqr(DO1)*sqr(DO4)*sqr(O2) +
                                               sqr(DO2)*(sqr(DO4)*(-16*sqr(DO1) +
                                               sqr(1 - 4*O1)) + 16*sqr(DO1)*sqr(O4)))));
          double rw1 = sqrt(abs(w1));
          dw1.first  *=0.5/rw1;
          dw1.second *=0.5/rw1;
          if (w1<0.) rw1*=-1.;
          double root2 = sqrt((4*sqr(DO2)*sqr(DO4)*sqr(O1) + sqr(DO1)*(sqr(DO4) *
                         sqr(1 - 2*O2) + 4*sqr(DO2)*sqr(O4)))/sqr(-1 + 4*O1 + 2*O2 + 2*O4));
          pair<double,double> dw2 = make_pair(DO1*DO2*DO4*(2*DO1*DO2*DO4 + (-1 + 4*O1 + 2*O2 + 2*O4)*root2) /
                                              ((-1 + 4*O1 + 2*O2 + 2*O4)*(4*sqr(DO2)*sqr(DO4)*sqr(O1) +
                                              sqr(DO1)*(sqr(DO4)*(-4*sqr(DO2) + sqr(1 - 2*O2)) + 4*sqr(DO2)*sqr(O4)))),
                                              DO1*DO2*DO4*(-2*DO1*DO2*DO4 + (-1 + 4*O1 + 2*O2 + 2*O4)*root2) /
                                              ((-1 + 4*O1 + 2*O2 + 2*O4)*(4*sqr(DO2)*sqr(DO4)*sqr(O1) +
                                              sqr(DO1)*(sqr(DO4)*(-4*sqr(DO2) +
                                              sqr(1 - 2*O2)) + 4*sqr(DO2)*sqr(O4)))));
          double rw2 = sqrt(abs(w2));
          dw2.first  *=0.5/rw2;
          dw2.second *=0.5/rw2;
          if (w2<0.) rw2*=-1.;
          double root3 = sqrt((4*sqr(DO2)*sqr(DO4)*sqr(O1) + sqr(DO1)*(4*sqr(DO4)*sqr(O2) +
                         sqr(DO2)*sqr(1 - 2*O4)))/sqr(-1 + 4*O1 + 2*O2 + 2*O4));
          pair<double,double> dw4 = make_pair(DO1*DO2*DO4*(2*DO1*DO2*DO4 + (-1 + 4*O1 + 2*O2 + 2*O4)*root3) /
                                              ((4*sqr(DO2)*sqr(DO4)*sqr(O1) + sqr(DO1)*(4*sqr(DO4)*sqr(O2) +
                                              sqr(DO2)*(-4*sqr(DO4) + sqr(1 - 2*O4)))) * (-1 + 4*O1 + 2*O2 + 2*O4)),
                                              DO1*DO2*DO4*(-2*DO1*DO2*DO4 + (-1 + 4*O1 + 2*O2 + 2*O4)*root3) /
                                              ((4*sqr(DO2)*sqr(DO4)*sqr(O1) + sqr(DO1)*(4*sqr(DO4)*sqr(O2) +
                                              sqr(DO2)*(-4*sqr(DO4) + sqr(1 - 2*O4)))) * (-1 + 4*O1 + 2*O2 + 2*O4)));
          double rw4 = sqrt(abs(w4));
          dw4.first  *=0.5/rw4;
          dw4.second *=0.5/rw4;
          if (w4<0.) rw4*=-1.;
          Estimate0DPtr h_w1;
          book(h_w1, 1, 3, 1);
          h_w1->set(rw1, dw1);
          Estimate0DPtr h_w2;
          book(h_w2, 1, 3, 2);
          h_w2->set(rw2, dw2);
          Estimate0DPtr h_w4;
          book(h_w4, 1, 3, 3);
          h_w4->set(rw4, dw4);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _n[3],_m[3][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2627838);

}
