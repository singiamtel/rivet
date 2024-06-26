// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > bbar and Upsilon(1,2,3) pi+pi-
  class BELLE_2015_I1336624 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1336624);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_hadrons, "TMP/c_hadrons", refData<YODA::BinnedEstimate<string>>(2,1,1));
      for(unsigned int ix=0;ix<3;++ix)
        book(_c_Ups[ix], "/TMP/c_Ups_"+toString(ix+1),refData<YODA::BinnedEstimate<string>>(1,1,1+ix));
      book(_c_muons[0], "/TMP/sigma_muons_1",refData<YODA::BinnedEstimate<string>>(2,1,1));
      book(_c_muons[1], "/TMP/sigma_muons_2",refData<YODA::BinnedEstimate<string>>(1,1,1));
      for(unsigned int ix=0;ix<2;++ix) {
        for (const string& en : _c_muons[ix].binning().edges<0>()) {
          const double end = std::stod(en)*MeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      // analyse the final state
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // intermediates
      bool isBottom(false);
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	// check for bottom hadrons
        if (PID::isBottomHadron(p.pid()) &&
            (p.abspid()%1000)/10!=55) {
	  isBottom = true;
	  break;
	}
	// upsilon + pi+pi-
	if(p.children().empty()) continue;
	if(p.pid() !=   553  &&
	   p.pid() != 100553 &&
	   p.pid() != 200553 ) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211) {
	    continue;
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(matched) {
	  if(nRes[211]==1 && nRes[-211]==1 ) {
	    if(p.pid()==553)
	      _c_Ups[0]->fill(_ecms[1]);
	    if(p.pid()==100553)
	      _c_Ups[1]->fill(_ecms[1]);
	    if(p.pid()==200553)
	      _c_Ups[2]->fill(_ecms[1]);
	  }
	}
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22]) {
	_c_muons[0]->fill(_ecms[0]);
	_c_muons[1]->fill(_ecms[1]);
      }
      // open bottom
      else if(isBottom) {
	_c_hadrons->fill(_ecms[0]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> tmp;
      book(tmp,2,1,1);
      divide(_c_hadrons,_c_muons[0],tmp);
      for(unsigned int ix=0;ix<3;++ix) {
        book(tmp,1,1,1+ix);
        divide(_c_Ups[ix],_c_muons[1],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons[2], _c_Ups[3];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1336624);

}
