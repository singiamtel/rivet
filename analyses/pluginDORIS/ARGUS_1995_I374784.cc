// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief tau -> e/mu nu nu
  class ARGUS_1995_I374784 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1995_I374784);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(ChargedFinalState(),"CFS");
      // hists
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double mtau=1.77686;
      Particles charged = apply<ChargedFinalState>(event,"CFS").particles();
      if (charged.size()!=4) vetoEvent;
      // find the charged lepton
      int il=-1;
      for (unsigned int ix=0; ix<charged.size(); ++ix) {
        if (charged[ix].abspid()==11 || charged[ix].abspid()==13) {
          if (il>0) vetoEvent;
          il=ix;
        }
      }
      if (il<0) vetoEvent;
      // cuts angle charged<90 each other and >90 lepton
      FourMomentum phad;
      Vector3 axisL = charged[il].mom().p3().unit();
      for (int ix=0; ix<4; ++ix) {
        if (ix==il) continue;
        phad += charged[ix];
        Vector3 axis = charged[ix].mom().p3().unit();
        if (axis.dot(axisL)>0) vetoEvent;
        for (int iy=ix+1; iy<4; ++iy) {
          if (iy==il) continue;
          if (axis.dot(charged[iy].mom().p3().unit())<0) vetoEvent;
        }
      }
      const double pp = sqrt(sqr(0.5*sqrtS())-sqr(mtau));
      Vector3 axis =-phad.p3().unit();
      FourMomentum ptau;
      ptau.setT(0.5*sqrtS());
      ptau.setX(pp*axis.x());
      ptau.setY(pp*axis.y());
      ptau.setZ(pp*axis.z());
      LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(ptau.betaVec());
      FourMomentum pl = boost.transform(charged[il].mom());
      if (charged[il].abspid()==11) {
      	_h[0]->fill(pl.p3().mod());
      }
      else {
      	_h[1]->fill(pl.p3().mod());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1995_I374784);

}
