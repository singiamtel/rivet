// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> D(*) D(*) K
  class BABAR_2008_I758167 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I758167);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<int>("PID", 9030443);
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 ||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }

    bool isK(int id) {
      return id==321 || id==311 || id==310 || id==130; 
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for(const Particle & B : ufs.particles()) {
	if(B.children().size()!=2) continue;
	int sign = B.pid()/B.abspid(), iD = B.abspid()-100;
	unsigned int imode=0;
	// boost to B rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B.momentum().betaVec());
	Particle resonance;
	if(B.children()[0].pid()== sign*10433 &&
	   B.children()[1].pid()==-sign*iD) {
	  imode=0;
	  resonance=B.children()[0];
	}
	else if(B.children()[1].pid()== sign*10433 &&
		B.children()[0].pid()==-sign*iD) {
	  imode=0;
	  resonance=B.children()[1];
	}
	else if(B.children()[0].pid()== 30443 &&
		isK(B.children()[1].pid())) {
	  imode=1;
	  resonance=B.children()[0];
	}
	else if(B.children()[1].pid()== 30443 &&
		isK(B.children()[0].pid())) {
	  imode=1;
	  resonance=B.children()[1];
	}
	else if(B.children()[0].pid()== _pid &&
		isK(B.children()[1].pid())) {
	  imode=2;
	  resonance=B.children()[0];
	}
	else if(B.children()[1].pid()== _pid &&
		isK(B.children()[0].pid())) {
	  imode=2;
	  resonance=B.children()[1];
	}
	else
	  continue;
	if(resonance.children().size()!=2) continue;
	FourMomentum pRes = boost1.transform(resonance.momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pRes.betaVec());
	Vector3 axis1 = -pRes.p3().unit();
	Particle child;
	if(imode==0) {
	  Particle Dstar;
	  if(isK(resonance.children()[1].pid()) &&
	     (resonance.children()[0].pid()==sign*413 || resonance.children()[0].pid()==sign*423))
	    child = resonance.children()[0];
	  else if(isK(resonance.children()[0].pid()) &&
		  (resonance.children()[1].pid()==sign*413 || resonance.children()[1].pid()==sign*423))
	    child = resonance.children()[1];
	  else
	    continue;
	}
	// B -> D psi(3770)
	else if(imode==1) {
	  if(resonance.children()[0].pid()!=-resonance.children()[1].pid()) continue;
	  if(resonance.children()[0].abspid()!=411 &&
	     resonance.children()[0].abspid()!=421) continue;
	  child = resonance.children()[0].pid()>0 ? resonance.children()[0] : resonance.children()[1];
	}
	else if(imode==2) {
	  if( (resonance.children()[0].abspid()==423 && resonance.children()[1].abspid()==421) ||
	      (resonance.children()[1].abspid()==423 && resonance.children()[0].abspid()==421)) {
	    child = resonance.children()[0].pid()>0 ? resonance.children()[0] : resonance.children()[1];
	  }
	  else continue;
	}
	Vector3 axis2 = boost2.transform(boost1.transform(child.momentum())).p3().unit();
	_h[imode]->fill(axis1.dot(axis2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I758167);

}
