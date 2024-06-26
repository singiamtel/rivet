// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > rho+ rho-
  class BABAR_2008_I789011 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I789011);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::abspid==213), "UFS");
      declare(FinalState(), "FS");
      // histos
      for(unsigned int ix=0;ix<5;++ix)
	book(_h[ix],1,1,1+ix);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // loop over rho mesons
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      Particle rhoP,rhoM;
      bool matched(false);
      for (const Particle& p :  ufs.particles()) {
	if(p.children().empty()) continue;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	matched=false;
	// check for antiparticle
	for (const Particle& p2 :  ufs.particles(Cuts::pid==-p.pid())) {
	  if(p2.children().empty()) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2==0) {
	    matched = true;
	    for(auto const & val : nRes2) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    // found rho+/-
	    if(matched) {
	      if(p.pid()>0) {
		rhoP = p;
		rhoM = p2;
	      }
	      else {
		rhoP = p2;
		rhoM = p;
	      }	
	      break;
	    }
	  }
	}
	if(matched) break;
      }
      if(!matched) vetoEvent;
      Particle piP,piM;
      matched = false;
      if (rhoP.children().size()!=2) vetoEvent;
      if (rhoM.children().size()!=2) vetoEvent;
      if(rhoP.children()[0].pid()==211 &&
	 rhoP.children()[1].pid()==111)
	piP = rhoP.children()[0];
      else if(rhoP.children()[1].pid()==211 &&
	      rhoP.children()[0].pid()==111)
	piP = rhoP.children()[1];
      else
	vetoEvent;
      if(rhoM.children()[0].pid()==-211 &&
	 rhoM.children()[1].pid()==111)
	piM = rhoM.children()[0];
      else if(rhoM.children()[1].pid()==-211 &&
	      rhoM.children()[0].pid()==111)
	piM = rhoM.children()[1];
      else
	vetoEvent;
      // boost to the rho+ rest frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(rhoP.momentum().betaVec());
      Vector3 e1z = rhoP.momentum().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      Vector3 axis1 = boost1.transform(piP.momentum()).p3().unit();
      double n1x(e1x.dot(axis1)),n1y(e1y.dot(axis1)),n1z(e1z.dot(axis1));
      _h[0]->fill(n1z);
      _h[2]->fill(atan2(n1y,n1x));
      // boost to the rho- frame
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(rhoM.momentum().betaVec());
      Vector3 axis2 = boost2.transform(piM.momentum()).p3().unit();
      double n2x(e1x.dot(axis2)),n2y(e1y.dot(axis2)),n2z(e1z.dot(axis2));
      _h[1]->fill(n2z);
      _h[3]->fill(atan2(n2y,n2x));
      _h[4]->fill(axis.dot(rhoP.momentum().p3().unit()));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I789011);

}
