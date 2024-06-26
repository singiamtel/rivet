#// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Charm cross sections 3.92 and 4.26 GeV
  class CLEOC_2008_I777917 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2008_I777917);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix)
        for(unsigned int iy=0;iy<3;++iy)
          book(_sigma_DD[ix][iy],1+ix,1,1+iy);
      for(unsigned int ix=0;ix<2;++ix) {
        book(_sigma_DDpi[ix],4,1,1+ix);
        book(_sigma_DDX [ix],5,1,1+ix);
      }
      book(_sigma_R[0],"TMP/hadron",refData<YODA::BinnedEstimate<int>>(6,1,1));
      book(_sigma_R[1],"TMP/muon",refData<YODA::BinnedEstimate<int>>(6,1,1));
      book(_sigma_cc,6,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22]) {
	_sigma_R[1]->fill(round(sqrtS()/MeV));
        return;
      }
      else
        _sigma_R[0]  ->fill(round(sqrtS()/MeV));
      // identified final state with D mesons
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
	int id1 = abs(p1.pid());
	if(id1 != 411 && id1 != 413 && id1 != 421 && id1 != 423 &&
	   id1 != 431 && id1 != 433)
	  continue;
	// check fs
	bool fs = true;
	for (const Particle & child : p1.children()) {
	  if(child.pid()==p1.pid()) {
	    fs = false;
	    break;
	  }
	}
	if(!fs) continue;
	// find the children
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	bool matched=false;
	int sign = p1.pid()/id1;
	// loop over the other fs particles
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
	  const Particle& p2 = ufs.particles()[iy];
	  fs = true;
	  for (const Particle & child : p2.children()) {
	    if(child.pid()==p2.pid()) {
	      fs = false;
	      break;
	    }
	  }
	  if(!fs) continue;
	  if(p2.pid()/abs(p2.pid())==sign) continue;
	  int id2 = abs(p2.pid());
	  if(id2 != 411 && id2 != 413 && id2 != 421 && id2 != 423 &&
	     id2 != 431 && id2 != 433)
	    continue;
	  if(!p2.parents().empty() && p2.parents()[0].pid()==p1.pid())
	    continue;
	  if((id1==411 || id1==421 || id1==431) && (id2==411 || id2==421 || id2==431 )) {
	    _sigma_DDX[1]->fill(round(sqrtS()/MeV));
            _sigma_cc->fill(round(sqrtS()/MeV));
          }
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2==0) {
	    matched=true;
	    for(auto const & val : nRes2) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      if(id1==411 && id2==411) {
		_sigma_DD[1][0]->fill(round(sqrtS()/MeV));
		_sigma_DDX[0]  ->fill(round(sqrtS()/MeV));
	      }
	      else if(id1==421&& id2==421) {
		_sigma_DD[0][0]->fill(round(sqrtS()/MeV));
		_sigma_DDX[0]  ->fill(round(sqrtS()/MeV));
	      }
	      else if(id1==431&& id2==431) {
		_sigma_DD[2][0]->fill(round(sqrtS()/MeV));
	      }
	      else if(id1==413 && id2==413) {
		_sigma_DD[1][2]->fill(round(sqrtS()/MeV));
	      }
	      else if(id1==423&& id2==423) {
		_sigma_DD[0][2]->fill(round(sqrtS()/MeV));
	      }
	      else if(id1==433&& id2==433) {
		_sigma_DD[2][2]->fill(round(sqrtS()/MeV));
	      }
	      else if((id1==421 && id2==423) ||
		      (id1==423 && id2==421)) {
		_sigma_DD[0][1]->fill(round(sqrtS()/MeV));
	      }
	      else if((id1==411 && id2==413) ||
		      (id1==413 && id2==411)) {
		_sigma_DD[1][1]->fill(round(sqrtS()/MeV));
	      }
	      else if((id1==431 && id2==433) ||
		      (id1==433 && id2==431)) {
		_sigma_DD[2][1]->fill(round(sqrtS()/MeV));
	      }
	    }
	  }
	  else if(ncount2==1) {
	    int ipi=0;
	    if(nRes2[111]==1 && nRes2[211]==0 && nRes[-211]==0 )
	      ipi = 111;
	    else if(nRes2[111]==0 && nRes2[211]==1 && nRes[-211]==0 )
	      ipi = 211;
	    else if(nRes2[111]==0 && nRes2[211]==0 && nRes[-211]==1 )
	      ipi =-211;
	    if(ipi==0) continue;
	    matched=true;
	    for(auto const & val : nRes2) {
	      if(val.first==ipi)
		continue;
	      else if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      bool Ddecay = false;
	      Particle mother = p1;
	      while (!mother.parents().empty()) {
		mother = mother.parents()[0];
		if(PID::isCharmMeson(mother.pid()) && mother.pid()!=p1.pid()) {
		  Ddecay = true;
		  break;
		}
	      }
	      mother = p2;
	      while (!mother.parents().empty()) {
		mother = mother.parents()[0];
		if(PID::isCharmMeson(mother.pid()) && mother.pid()!=p1.pid()) {
		  Ddecay = true;
		  break;
		}
	      }
	      if(Ddecay) continue;
	      if((id1==413 || id1==423 ) &&
		 (id2==413 || id2==423 )) {
		_sigma_DDpi[1]->fill(round(sqrtS()/MeV));
	      }
	      else if((id1==411 || id1==421 ) &&
		      (id2==413 || id2==423 )) {
		_sigma_DDpi[0]->fill(round(sqrtS()/MeV));
	      }
	      else if((id1==413 || id1==423 ) &&
		      (id2==411 || id2==421 )) {
		_sigma_DDpi[0]->fill(round(sqrtS()/MeV));
	      }
	    }
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // cross sections
      double fact = crossSection()/picobarn/ sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix)
        for(unsigned int iy=0;iy<3;++iy)
          scale(_sigma_DD[ix][iy],fact);
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_sigma_DDpi[ix],fact);
        scale(_sigma_DDX [ix],fact*1e-3); // this one in nb
      }
      BinnedEstimatePtr<int> tmp;
      book(tmp,6,1,2);
      divide(_sigma_R[0],_sigma_R[1],tmp);
      scale(_sigma_cc,fact*1e-3);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _sigma_DD[3][3];
    BinnedHistoPtr<int> _sigma_DDpi[2];
    BinnedHistoPtr<int> _sigma_DDX[2];
    BinnedHistoPtr<int> _sigma_R[2],_sigma_cc;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2008_I777917);


}
