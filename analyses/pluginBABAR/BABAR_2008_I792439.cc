// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> phi K pi
  class BABAR_2008_I792439 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I792439);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(333);
      B0.addStable(310);
      B0.addStable(111);
      declare(B0, "B0");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
	book(_p[ix][0],1,2,2+ix);
	book(_p[ix][1],"TMP/norm_"+toString(ix));
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_mass[ix][iy],2,1+ix,1+iy);
	  book(_h_angle[ix][iy],3+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1},{-211,1}, { 333,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1},{ 211,1}, { 333,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 310,1},{ 111,1}, { 333,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1,imode=-1;
      	if (B0.modeMatches(ix,3,mode1)) {
	  imode = 0;
      	  sign  = 1;
	}
	else if(B0.modeMatches(ix,3,mode1CC)) {
	  imode = 0;
      	  sign  =-1;
	}
	else if(B0.modeMatches(ix,3,mode2)) {
	  imode = 1;
      	  sign  = 1;
	}
      	else
      	  continue;
	// particles
      	const Particle & KK  = B0.decayProducts()[ix].at(imode==0 ?  321*sign : 310)[0];
      	const Particle & pi  = B0.decayProducts()[ix].at(imode==0 ? -211*sign : 111)[0];
	const Particle & phi = B0.decayProducts()[ix].at( 333     )[0];
	// children of the phi
	if(phi.children().size()!=2) continue;
	if(phi.children()[0].abspid()!=321) continue;
	if(phi.children()[0].pid()!=-phi.children()[1].pid()) continue;
	double mKpi = (KK.momentum()+pi.momentum()).mass();
	_h_mass[imode][0]->fill(mKpi);
	_h_mass[imode][1]->fill(phi.mass());
	Particle Kp1 = phi.children()[0];
	Particle Km1 = phi.children()[1];
	// B0 frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum pKstar = boost1.transform(KK.momentum()+pi.momentum());
	FourMomentum pPhi   = boost1.transform(phi.momentum());
	// K pi helicity angle
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	FourMomentum pKp = boost2.transform(boost1.transform(KK.momentum()));
	Vector3 axis1 = pKstar.p3().unit();
	double cTheta1 = axis1.dot(pKp.p3().unit());
	// phi helicity angle
	LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
	FourMomentum pKp1 = boost3.transform(boost1.transform(Kp1.momentum()));
	Vector3 axis2 = pPhi.p3().unit();
	double cTheta2 = axis2.dot(pKp1.p3().unit());
	if(mKpi>.75 && mKpi<1.05) {
	  _h_angle[0][0]->fill(cTheta1);
	  _h_angle[0][1]->fill(cTheta2);
	  _p[0][0]->fill(-0.5*(1-5.*sqr(cTheta2)));
	  _p[0][1]->fill();
	}
	else if(mKpi>1.13 && mKpi<1.53) {
	  _h_angle[1][0]->fill(cTheta1);
	  _h_angle[1][1]->fill(cTheta2);
	  _p[1][0]->fill(-0.5*(1-5.*sqr(cTheta2)));
	  _p[1][1]->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_p[ix][0], 1./ *_p[ix][1]);
	for(unsigned int iy=0;iy<2;++iy) {
	  normalize(_h_mass[ix][iy] ,1.,false);
	  normalize(_h_angle[ix][iy],1.,false);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2][2],_h_angle[2][2];
    CounterPtr _p[2][2];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I792439);

}
