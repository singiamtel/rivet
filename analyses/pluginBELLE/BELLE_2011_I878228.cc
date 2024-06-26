// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > Ds(*)+ Ds*()-
  class BELLE_2011_I878228 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2011_I878228);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(  _c["DpDm"], "/TMP/sigma_DpDm"  , refData(1,1,1));
      book( _c["DpDmS"], "/TMP/sigma_DpDmS" , refData(1,1,2));
      book(_c["DpSDmS"], "/TMP/sigma_DpSDmS", refData(1,1,3));
      book(   _c["All"], "/TMP/sigma_All"   , refData(2,1,1));
      book(    _c["mu"], "/TMP/sigma_mu"    , refData(2,1,1));
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
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
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        _c["mu"]->fill(sqrtS());
        return;
      }
      // unstable charm analysis

      Particles ds = apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==431 || Cuts::abspid==433);
      for (unsigned int ix=0; ix<ds.size(); ++ix) {
       	const Particle& p1 = ds[ix];
        int id1 = abs(p1.pid());
      	// check fs
      	bool fs = true;
      	for (const Particle & child : p1.children()) {
      	  if (child.pid()==p1.pid()) {
      	    fs = false;
      	    break;
      	  }
      	}
      	if (!fs) continue;
      	// find the children
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
      	bool matched=false;
       	int sign = p1.pid()/id1;
        // loop over the other fs particles
        for (unsigned int iy=ix+1; iy<ds.size(); ++iy) {
          const Particle& p2 = ds[iy];
          fs = true;
          for (const Particle & child : p2.children()) {
            if (child.pid()==p2.pid()) {
              fs = false;
              break;
            }
          }
          if (!fs) continue;
          if (p2.pid()/abs(p2.pid())==sign) continue;
          int id2 = abs(p2.pid());
          if (!p2.parents().empty() && p2.parents()[0].pid()==p1.pid()) {
            continue;
          }
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2!=0) continue;
          matched=true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            if (id1==431 && id2==431) {
              _c["DpDm"]->fill(sqrtS());
            }
            else if (id1==433 && id2==433) {
              _c["DpSDmS"]->fill(sqrtS());
            }
            else if ((id1==431 && id2==433) || (id1==433 && id2==431)) {
              _c["DpDmS"]->fill(sqrtS());
            }
            _c["All"]->fill(sqrtS());
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c, crossSection()/ sumOfWeights()/nanobarn);
      size_t idx = 0;
      for (const string& s : { "DpDm"s, "DpDmS"s, "DpSDmS"s}) {
        Estimate1DPtr mult;
        book(mult, 1, 1, ++idx);
        barchart(_c[s],mult);
      }
      // R
      Estimate1DPtr mult;
      book(mult, 2,1,1);
      divide(_c["All"s],_c["mu"s],mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2011_I878228);

}
