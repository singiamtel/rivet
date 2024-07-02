// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  B -> J/psi / psi(2S) K pi
  class BABAR_2009_I801589 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I801589);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 ||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      BB.addStable(443);
      BB.addStable(100443);
      declare(BB, "BB");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_Kpi2[ix],2,1,1+ix);
	book(_h_angle[ix],4,1,1+ix);
	book(_c[ix],"TMP/c_"+toString(ix));
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_Kpi[ix][iy],1,1+ix,1+iy);
	}
        book(_b[ix],{0.,.795,.995,1.332,1.532,10.});
	for(unsigned int iy=0;iy<5;++iy)
	  book(_b[ix]->bin(iy+1),3,1+ix,1+iy);
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_h_PsiPi[ix][iy],5,1+ix,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1},{-211,1}, {    443,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1},{ 211,1}, {    443,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 321,1},{-211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-321,1},{ 211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 310,1},{-211,1}, {    443,1}};
      static const map<PdgId,unsigned int> & mode3CC = { { 310,1},{ 211,1}, {    443,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 310,1},{-211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode4CC = { { 310,1},{ 211,1}, { 100443,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int sign = 1,iK(0),iPsi(0);
	if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
	  sign=1; iK = 321; iPsi=443;
	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
       	  sign=-1; iK=-321; iPsi=443;
       	}
	else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	  sign=1; iK = 321; iPsi=100443;
	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
       	  sign=-1; iK=-321; iPsi=100443;
       	}
	else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3)) {
	  sign=1; iK = 310; iPsi=443;
	}
      	else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3CC)) {
       	  sign=-1; iK= 310; iPsi=443;
	}
	else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode4)) {
	  sign=1; iK = 310; iPsi=100443;
	}
      	else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode4CC)) {
       	  sign=-1; iK= 310; iPsi=100443;
       	}
      	else
      	  continue;
	_c[iPsi/100000]->fill();
       	const Particle & Kp  = BB.decayProducts()[ix].at( iK      )[0];
       	const Particle & pim = BB.decayProducts()[ix].at(-211*sign)[0];
       	const Particle & psi = BB.decayProducts()[ix].at( iPsi    )[0];
	FourMomentum pKpi = Kp.momentum()+pim.momentum(); 
       	double mKpi  = pKpi.mass();
	_h_Kpi[BB.decaying()[ix].abspid()%100/10-1][iPsi/100000]->fill(mKpi);
	_h_Kpi2[iPsi/100000]->fill(mKpi);
	double mPsiPi = (psi.momentum()+pim.momentum()).mass();
	_b[iPsi/100000]->fill(mKpi,mPsiPi);
	_h_PsiPi[iPsi/100000][0]->fill(mPsiPi);
	// helicity angle
	if(mKpi<0.795 || (mKpi>0.995&&mKpi<1.332) || mKpi>1.532) {
	  _h_PsiPi[iPsi/100000][2]->fill(mPsiPi);
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].momentum().betaVec());
	  pKpi = boost1.transform(pKpi);
	  Vector3 axis1 = pKpi.p3().unit();
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKpi.betaVec());
	  FourMomentum ppi = boost3.transform(boost1.transform(pim.momentum()));
	  double cPi = ppi.p3().unit().dot(axis1);
	  _h_angle[iPsi/100000]->fill(cPi);
	}
	else {
	  _h_PsiPi[iPsi/100000][1]->fill(mPsiPi);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_Kpi2[ix],1.,false);
	normalize(_h_angle[ix],1.,false);
	for(unsigned int iy=0;iy<2;++iy) {
	  normalize(_h_Kpi[ix][iy],1.,false);
	}
	for(unsigned int iy=0;iy<3;++iy) {
	  scale(_h_PsiPi[ix][iy],1./ *_c[ix]);
	}
        scale(_b[ix], 1./ *_c[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi[2][2],_h_Kpi2[2],_h_angle[2],_h_PsiPi[2][3];
    Histo1DGroupPtr _b[2];
    CounterPtr _c[2];
    /// @}
  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I801589);

}
