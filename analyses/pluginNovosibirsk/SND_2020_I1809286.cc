// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > pi+pi-pi0
  class SND_2020_I1809286 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2020_I1809286);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
        book(_sigma_total[ix],1+ix,1,1);
        for (const string& en : _sigma_total[ix].binning().edges<0>()) {
          double end = std::stod(en)*GeV;
          if(isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      for(unsigned int ix=0;ix<3;++ix)
        book(_sigma_res[ix], "/TMP/c_res_"+toString(ix),refData(3,1,1+ix));
      
      if(inRange(sqrtS()/GeV,1.42,1.48)) {
	book(_h_x,4,1,1);
	book(_h_m,4,1,3);
      }
      else if(inRange(sqrtS()/GeV,1.65,1.68)) {
	book(_h_x,4,1,2);
	book(_h_m,4,1,4);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
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

      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      Particle pip,pim;
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
	if(p.pid()     == 211) pip=p;
	else if(p.pid()==-211) pim=p;
      }
      if(ntotal!=3) vetoEvent;
      if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==1) {
	_sigma_total[0]->fill(_ecms[0]);
	_sigma_total[1]->fill(_ecms[1]);
      }
      else
	vetoEvent;
      if(_h_x) {
	_h_x->fill(pip.momentum().p()/sqrtS());
	_h_x->fill(pim.momentum().p()/sqrtS());
	_h_m->fill((pip.momentum()+pim.momentum()).mass()/MeV);
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==223 or
					     Cuts::abspid==113    or Cuts::abspid==213 or
					     Cuts::abspid==100113 or Cuts::abspid==100213)) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=1) continue;
	int idOther = 211;
	if(p.pid()==223 || p.pid()==113 || p.pid()==100113)
	  idOther = 111;
	else if(p.pid()==213 || p.pid()==100213)
	  idOther = -211;
	bool matched=true;
	for(auto const & val : nRes) {
	  if(val.first==idOther ) {
	    if(val.second !=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(!matched) continue;
	if(matched) {
	  if(p.pid()==223)
	    _sigma_res[2]->fill(sqrtS());
	  else if(p.pid()==213 || p.pid()==113)
	    _sigma_res[0]->fill(sqrtS());
	  else
	    _sigma_res[1]->fill(sqrtS());
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_h_x) {
        normalize(_h_x,1.,false);
        normalize(_h_m,1.,false);
      }
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix)
        scale(_sigma_total[ix],fact);
      for(unsigned int ix=0;ix<3;++ix) {
        scale(_sigma_res[ix],fact);
        Estimate1DPtr tmp;
        book(tmp,3,1,1+ix);
        barchart(_sigma_res[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma_total[2];
    string _ecms[2];
    Histo1DPtr _sigma_res[3];
    Histo1DPtr _h_x,_h_m;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2020_I1809286);

}
