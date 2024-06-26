// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > J/psi X(3940)
  class BELLE_2007_I686580 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I686580);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<int>("PID", 9010441);
      // projections
      declare("FS",FinalState());
      declare("UFS",UnstableParticles(Cuts::pid==443 || Cuts::pid==_pid));
      // histograms
      book(_h,1,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount,
		      unsigned int & nCharged) {
      for( const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
          if(PID::isCharged(p.pid())) ++nCharged;
        }
        else {
          findChildren(child,nRes,ncount,nCharged);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // loop over J/psi
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      bool matched=false;
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        unsigned int nCharged=0;
        findChildren(p,nRes,ncount,nCharged);
        // X(3940)
      	for (const Particle& p2 : ufs.particles(Cuts::pid==_pid)) {
      	  map<long,int> nResB = nRes;
      	  int ncountB = ncount;
      	  unsigned int nChargedB=0;
      	  findChildren(p2,nResB,ncountB,nChargedB);
      	  if (ncountB!=0) continue;
      	  matched = true;
      	  for (const auto& val : nResB) {
      	    if (val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  if (matched) {
            if (nCharged>2) _h->fill("10.6"s);
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection()/ sumOfWeights() /femtobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    BinnedHistoPtr<string> _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I686580);

}
