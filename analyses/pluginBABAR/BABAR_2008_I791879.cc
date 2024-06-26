// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> Lambdac+ pbar pi-
  class BABAR_2008_I791879 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I791879);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==PID::BPLUS);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 4122);
      BP.addStable(-4122);
      declare(BP, "BP");
      // histos
      book(_h[0],"TMP/hA",refData(1,1,1));
      book(_h[1],2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-4122,1}, { 2212,1}, { 211,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 4122,1}, {-2212,1}, {-211,1} };
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	// select right decay mode
	int sign = 1;
	if      (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode  )) sign = 1;
	else if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) sign =-1;
	else continue;
	const Particle & lam = BP.decayProducts()[ix].at(-4122*sign)[0];
	const Particle & pp  = BP.decayProducts()[ix].at( 2212*sign)[0];
	const Particle & pip = BP.decayProducts()[ix].at(  211*sign)[0];
	_h[0]->fill((lam.momentum()+pp.momentum()).mass());
	if(BP.decaying()[ix].children().size()==2 &&
	   (BP.decaying()[ix].children()[0].abspid()==4112 ||
	    BP.decaying()[ix].children()[1].abspid()==4112)) {
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].momentum().betaVec());
	  FourMomentum pSigma  = boost1.transform(lam.momentum()+pip.momentum());
	  Vector3 axis = -pSigma.p3().unit();
	  FourMomentum pLambda = boost1.transform(lam.momentum());
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pSigma.betaVec());
	  _h[1]->fill(axis.dot(boost2.transform(pLambda).p3().unit()));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // phase space volumes in the different bins (PR has python program which computed these)
      vector<double> phsp = {0.1351399774095767,0.24107906667337267,0.3088204600383566,0.36206863378662535,
			     0.40659958403332336,0.4450070368053994,0.47872719863744506,0.5086545750109557,
			     0.5353906093802568,0.5593621540798246,0.5808845790744056,0.6001982286826937,
			     0.6174908272498764,0.6329119405055217,0.6465826880700479,0.6586024850910634,
			     0.6690538534165392,0.67800593732165,0.685517125504117,0.6916370414215524,
			     0.6964080775759626,0.6998665941923083,0.7020438666137178,0.7029668415116569,
			     0.7026587454114315,0.701139577428797,0.6984265098419529,0.6945342141075089,
			     0.6894751254590167,0.6832596558172149,0.6758963620528058,0.6673920744258324,
			     0.657751988084915,0.6469797186756903,0.6350773212151493,0.6220452692510539,
			     0.6078823887065968,0.5925857373802724,0.5761504163391165,0.5585692926468113};
      // normalize the histograms
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
      // convert first plot to scatter and normalize to phase space volume in bin
      // convert to scatter
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_h[0],tmp);
      double step = 0.0344;
      // divide by phase space volume
      for(unsigned int ip=0;ip<tmp->numBins();++ip)
	tmp->bin(ip+1).scale(1./phsp[ip]/step);
    }
    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I791879);

}
