// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief mu+mu- and tau+tau- in gamma gamma
  class L3_2004_I645127 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_2004_I645127);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // get the mode and options
      _mode =0;
      if( getOption("PROCESS") == "EE" ) _mode = 0;
      else if( getOption("PROCESS") == "GG") _mode = 1;

      // Initialise and register projections
      if(_mode==0) {
	declare(GammaGammaKinematics(), "Kinematics");
        declare(GammaGammaFinalState(), "FS");
        declare(UnstableParticles(),"UFS");
        // Book histos
        book(_c_sigma_mu1, 1,1,1);
        book(_c_sigma_mu2, 1,1,2);
        book(_c_sigma_tau, 2,1,1);
      }
      else if(_mode==1) {
        declare(FinalState(), "FS");
        book(_sigma,"TMP/sigma",refData(3,1,1));
      }

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
    	if(child.children().empty()) {
    	  --nRes[child.pid()];
    	  --ncount;
    	}
    	else
    	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // stuff for e+e- collisions
      double W2 = sqr(sqrtS());
      if(_mode==0) {
	const GammaGammaKinematics& kin = apply<GammaGammaKinematics>(event, "Kinematics");
	W2 = kin.W2();
	if(W2<9.*sqr(GeV) ) vetoEvent;
      }
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      bool fiducal = true;
      for(const Particle & p : fs.particles()) {
     	nCount[p.pid()] += 1;
     	++ntotal;
	if(abs(p.pid())==13) {
	  if(abs(cos(p.momentum().polarAngle()))>0.8) fiducal = false;
	}
      }
      if( nCount[-13]==1 && nCount[13]==1 && ntotal==2+nCount[22]) {
	if(W2<1600.*sqr(GeV) && _c_sigma_mu1) {
          _c_sigma_mu2->fill(round(sqrtS()));
          if(fiducal) _c_sigma_mu1->fill(round(sqrtS()));
        }
        if(_sigma) _sigma->fill(sqrtS());
      }
      if(_mode==1) return;
      bool foundTauPlus = false, foundTauMinus = true;
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
    	// find the taus
     	if(abs(p.pid())==15) {
	  if(p.pid()== 15) foundTauMinus=true;
	  if(p.pid()==-15) foundTauPlus =true;
     	  findChildren(p,nCount,ntotal);
	}
      }
      if(!foundTauPlus || !foundTauMinus)
	vetoEvent;
      bool matched = true;
      for(auto const & val : nCount) {
	if(val.first==22) {
	  continue;
	}
	else if(val.second!=0) {
	  matched = false;
	  break;
	}
      }
      if(matched)
	_c_sigma_tau->fill(round(sqrtS()));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // prefactor for the cross sections
      double fact  = crossSection()/picobarn/sumOfWeights();
      if(_mode==0) {
        scale(_c_sigma_mu1,fact);
        scale(_c_sigma_mu2,fact);
        scale(_c_sigma_tau,fact);
      }
      else {
        fact /= 1000.;
        scale(_sigma,fact);
        for(unsigned int iy=1;iy<6;++iy) {
          Estimate1DPtr tmp;
          book(tmp,3,1,iy);
          barchart(_sigma,tmp);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _c_sigma_mu1,_c_sigma_mu2,_c_sigma_tau;
    Histo1DPtr _sigma;
    unsigned int _mode;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(L3_2004_I645127);


}
