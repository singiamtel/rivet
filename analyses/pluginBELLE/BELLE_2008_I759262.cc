// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > J/psi X(3940) or X(4160)
  class BELLE_2008_I759262 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I759262);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid[0] = getOption<int>("PID3940", 9010441);
      _pid[1] = getOption<int>("PID4160", 9010445);
      // projections
      declare("FS",FinalState());
      declare("UFS",UnstableParticles(Cuts::pid==443 ||
                                      Cuts::pid==_pid[0] ||
                                      Cuts::pid==_pid[1]));
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount,
		      unsigned int & nCharged) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
          if (PID::isCharged(p.pid())) ++nCharged;
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
      	// X(3940)/X(4160)
       	for (const Particle& p2 : ufs.particles(Cuts::pid!=443)) {
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
            if (p2.children().size()!=2) break;
            if (p2.children()[0].pid()*p2.children()[1].pid()>0) break;
            if (p2.pid()==_pid[0]) {
              if ((p2.children()[0].abspid()==413 && p2.children()[1].abspid()==411) ||
                  (p2.children()[1].abspid()==413 && p2.children()[0].abspid()==411) ||
                  (p2.children()[0].abspid()==423 && p2.children()[1].abspid()==421) ||
                  (p2.children()[1].abspid()==423 && p2.children()[0].abspid()==421))
                _h[0]->fill("10.6"s);
              }
              else {
              if (p2.children()[0].pid()==-p2.children()[1].pid() &&
                 (p2.children()[0].abspid()==413 || p2.children()[0].abspid()==423))
                _h[1]->fill("10.6"s);
              }
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
    int _pid[2];
    BinnedHistoPtr<string> _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2008_I759262);

}
