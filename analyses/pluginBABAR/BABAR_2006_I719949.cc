// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e -> rho 0 rho and rho0 phi
  class BABAR_2006_I719949 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I719949);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::pid==113 or
				Cuts::pid==333), "UFS");
      declare(FinalState(), "FS");
      // histos
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_hel  [ix],3,1,1+ix);
	if(ix==2) continue;
	book(_h_sigma[ix],1,1,1+ix);
	book(_h_prod [ix],2,1,1+ix);
      }
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
      const Particles vMesons = apply<UnstableParticles>(event, "UFS").particles();
      Particle vectors[2];
      bool matched(false);
      for (unsigned int ix=0;ix<vMesons.size();++ix) {
	if(vMesons[ix].children().empty()) continue;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(vMesons[ix],nRes,ncount);
	matched=false;
	for (unsigned int iy=ix+1;iy<vMesons.size();++iy) {
	  if(vMesons[iy].children().empty()) continue;
	  if(vMesons[ix].pid()==333 && vMesons[iy].pid()==333) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(vMesons[iy],nRes2,ncount2);
	  if(ncount2==0) {
	    matched = true;
	    for(auto const & val : nRes2) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      vectors[0] = vMesons[ix];
	      vectors[1] = vMesons[iy];
	      break;
	    }
	  }
	}
	if(matched) break;
      }
      if(!matched) vetoEvent;
      if(vectors[0].pid()==333) swap(vectors[0],vectors[1]);
      if (vectors[0].children().size()!=2) vetoEvent;
      if (vectors[1].children().size()!=2) vetoEvent;
      double cTheta = abs(axis.dot(vectors[0].momentum().p3().unit()));
      if(cTheta>0.8) vetoEvent;
      if(vectors[0].pid()==vectors[1].pid()) {
	_h_sigma[0]->fill("10.58"s);
	_h_prod [0]->fill(cTheta);
      }
      else {
	_h_sigma[1]->fill("10.58"s);
	_h_prod [1]->fill(cTheta);
      }
      // helicity angles
      double cHel[2];
      for(unsigned int ix=0;ix<2;++ix) {
	int iMeson = vectors[ix].pid() == 113 ? 211 : 321;
	Particle mP;
	if(vectors[ix].children()[0].pid()== iMeson &&
	   vectors[ix].children()[1].pid()==-iMeson)
	  mP = vectors[ix].children()[0];
	else if(vectors[ix].children()[1].pid()== iMeson &&
		vectors[ix].children()[0].pid()==-iMeson)
	  mP = vectors[ix].children()[1];
	else
	  vetoEvent;
	// boost to the rho+ rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(vectors[ix].momentum().betaVec());
	Vector3 e1z = vectors[ix].momentum().p3().unit();
	Vector3 axis1 = boost1.transform(mP.momentum()).p3().unit();
	cHel[ix] = e1z.dot(axis1);
      }
      if(vectors[0].pid()==vectors[1].pid()) {
	_h_hel[0]->fill(cHel[0]);
	_h_hel[0]->fill(cHel[1]);
      }
      else {
	_h_hel[1]->fill(cHel[1]);
	_h_hel[2]->fill(cHel[0]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/sumOfWeights()/femtobarn;
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h_hel  [ix],1.,false);
	if(ix==2) continue;
	scale(_h_sigma[ix],fact);
	normalize(_h_prod [ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_sigma[2];
    Histo1DPtr _h_prod[2],_h_hel[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I719949);

}
