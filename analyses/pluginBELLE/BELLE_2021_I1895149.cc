// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> X_u l nu
  class BELLE_2021_I1895149 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1895149);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histograms
      for(unsigned int ix=0;ix<6;++ix) {
	book(_h_direct[ix],1+ix,1,1);
	book(_h_forward[ix],"TMP/h_"+toString(ix+1),refData(7+ix,1,1));
      }
      book(_nB,"/TMP/nB");
    }

    void findDecayProducts(Particle parent, Particles & em, Particles & ep,
			   Particles & nue, Particles & nueBar, bool & charm) {
      for(const Particle & p : parent.children()) {
	if(PID::isCharmHadron(p.pid())) {
	  charm=true;
	}
	else if(p.pid() == PID::EMINUS) {
	  em.push_back(p);
	}
	else if(p.pid() == PID::EPLUS) {
	  ep.push_back(p);
	}
	else if(p.pid() == PID::NU_E || p.pid()==PID::NU_MU) {
	  nue.push_back(p);
	}
	else if(p.pid() == PID::NU_EBAR || p.pid()==PID::NU_MUBAR) {
	  nueBar.push_back(p);
	}
	else if(PID::isBottomHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar,charm);
	}
	else if(!PID::isHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar,charm);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find and loop over Upslion(4S)
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
      	for(const Particle & p2 : p.children()) {
      	  if(p2.abspid()!=511 && p2.abspid()!=521) continue;
	  _nB->fill();
	  bool charm = false;
	  Particles em,ep,nue,nueBar;
	  findDecayProducts(p2,em,ep,nue,nueBar,charm);
	  if(charm) continue;
	  FourMomentum pl,pnu;
	  if(em.size()==1 && nueBar.size()==1) {
	    pl  = em[0].momentum();
	    pnu = nueBar[0].momentum();
	  }
	  else if(ep.size()==1 && nue.size()==1) {
	    pl  = ep[0].momentum();
	    pnu = nue[0].momentum();
	  }
	  else
	    continue;
      	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p2.momentum().betaVec());
	  pl  = boost.transform(pl );
	  pnu = boost.transform(pnu);
	  FourMomentum pB = boost.transform(p2.momentum());
	  FourMomentum q = pl+pnu;
	  FourMomentum pX = pB-q;
	  double p3 = pX.p();
	  _h_forward[0]->fill(pl.E());
	  _h_forward[1]->fill(q.mass2());
	  _h_forward[2]->fill(pX.mass());
	  _h_forward[3]->fill(pX.mass2());
	  _h_forward[4]->fill(pX.E()-p3);
	  _h_forward[5]->fill(pX.E()+p3);
	  if(pl.E()>1) {
	    _h_direct[0]->fill(pl.E());
	    _h_direct[1]->fill(q.mass2());
	    _h_direct[2]->fill(pX.mass());
	    _h_direct[3]->fill(pX.mass2());
	    _h_direct[4]->fill(pX.E()-p3);
	    _h_direct[5]->fill(pX.E()+p3);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<6;++ix) {
        // unfolded dist, scale by 1/2 /no of B's (2 as using e and mu modes)
        scale(_h_direct[ix], 0.5/ *_nB);
        // forward folding scale to BELLE no of B's
        scale(_h_forward[ix], 2.*771.58e6/ *_nB);
        // get the efficiency product and divide by it
        unsigned int iloc = ix<2 ? 3+ix : (ix<4 ? ix-1 : ix+1);
        Estimate1D eff    = refData<YODA::Estimate1D>(iloc+24,1,1);
        Estimate2D matrix = refData<YODA::Estimate2D>(  19+ix,1,1);
        // scatter for the result
        Estimate1DPtr corrected;
        book(corrected,ix+7,1,1);
        vector<double> val(_h_forward[ix]->numBins(),0.),err(_h_forward[ix]->numBins(),0.);
        // first divide by eff
        for (unsigned int iy=0;iy<_h_forward[ix]->numBins();++iy) {
          val[iy] = _h_forward[ix]->bin(iy+1).sumW()/eff.bin(iy+1).val();
          double relE = eff.bin(iy+1).totalErrAvg()/eff.bin(iy+1).val();
          err[iy] =val[iy]*sqrt(sqr(relE) + sqr(_h_forward[ix]->bin(iy+1).relErrW()));
        }
        vector<double> val2(_h_forward[ix]->numBins(),0.),err2(_h_forward[ix]->numBins(),0.);
        for (unsigned int iy=0;iy<_h_forward[ix]->numBins();++iy) {
          for (unsigned int iz=0;iz<_h_forward[ix]->numBins();++iz) {
            double corr  = matrix.bin((_h_forward[ix]->numBins()+2)*(iz+1)+iy+1).val()/100.;
            double ecorr = matrix.bin((_h_forward[ix]->numBins()+2)*(iz+1)+iy+1).totalErrAvg()/100.;
            val2[iy] += corr*val[iz];
            err2[iy] += sqr(ecorr*val[iz]) + sqr(corr*err[iz]);
          }
          err2[iy]  = val2[iy]*sqrt(err2[iy]/sqr(val2[iy]) +  sqr(9.78/771.58));
        }
        for (unsigned int ibin=0; ibin<_h_forward[ix]->numBins(); ++ibin) {
          const double dy = sqrt(err[ibin]);
          corrected->bin(ibin+1).set(val[ibin], dy);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_direct [6];
    Histo1DPtr _h_forward[6];
    CounterPtr _nB;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1895149);

}
