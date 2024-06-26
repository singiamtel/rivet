// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B0 -> D*+ D*-
  class BABAR_2005_I686355 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I686355);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      // histograms
      book(_p[0],1,1,1);
      book(_p[1],"TMP/wgt");
      book(_h,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles B0 = apply<UnstableParticles>(event, "UFS").particles();
      for(const Particle & p : B0) {
	if(p.children().size()!=2) continue;
	if(p.children()[0].pid()!=-p.children()[1].pid()) continue;
	if(p.children()[0].abspid()!=413) continue;
	Particle Dp = p.children()[0];
	Particle Dm = p.children()[1];
	if     (p.pid()>0 && Dp.pid()<0) swap(Dp,Dm);
	else if(p.pid()<0 && Dp.pid()>0) swap(Dp,Dm);
	// find children of the D mesons
	Particle pip,pim;
	if(Dp.children().size()!=2) continue;
	if(Dp.children()[0].abspid()==PID::PIPLUS &&
	   Dp.children()[1].abspid()==PID::D0)
	  pip = Dp.children()[0];
	else if (Dp.children()[1].abspid()==PID::PIPLUS &&
		 Dp.children()[0].abspid()==PID::D0)
	  pip = Dp.children()[1];
	else
	  continue;
	if(Dm.children().size()!=2) continue;
	if(Dm.children()[0].abspid()==PID::PIPLUS &&
	   Dm.children()[1].abspid()==PID::D0) {
	  pim = Dm.children()[0];
	}
	else if (Dm.children()[1].abspid()==PID::PIPLUS &&
		 Dm.children()[0].abspid()==PID::D0) {
	  pim = Dm.children()[1];
	}
	else
	  continue;
	// boost to rest frame
	LorentzTransform boostB  = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pDstarPlus  = boostB.transform(Dp.momentum());
	FourMomentum pDstarMinus = boostB.transform(Dm.momentum());
	LorentzTransform boostDp = LorentzTransform::mkFrameTransformFromBeta(pDstarPlus .betaVec());
	FourMomentum ppip = boostDp.transform(boostB.transform(pip.momentum()));
	LorentzTransform boostDm = LorentzTransform::mkFrameTransformFromBeta(pDstarMinus.betaVec());
	FourMomentum ppim = boostDm.transform(boostB.transform(pim.momentum()));
	Vector3 axisX = pDstarPlus.p3().unit();
	// y and z axis
	ppim = boostDp.transform(boostB.transform(pim.momentum()));
	Vector3 axisZ = axisX.cross(ppim.p3()).unit();
	// cThetaTr
	double cThetaTr = axisZ.dot(ppip.p3().unit());
	_h->fill(cThetaTr);
	_p[0]->fill(0.5*(5.*sqr(cThetaTr)-1));
        _p[1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
      scale(_p[0], 1./ *_p[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _p[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2005_I686355);

}
