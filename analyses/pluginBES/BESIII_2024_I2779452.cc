// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for  $e^+e^-\to\omega \chi_{c(1,2)}$ and $X(3872)$ for $\sqrt{s}=4.66$ to 4.95 GeV
  class BESIII_2024_I2779452 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2779452);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
        book(_h[ix],1+ix,1,1);
      for(const string& en : _h[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren2(const Particle& p, map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.pid()==443 ) {
          nRes[child.pid()] += 1;
          ++ncount;
        }
        else if (child.children().empty()) {
          ++nRes[child.pid()];
          ++ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    void findChildren(const Particle& p, map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
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
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p: fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& XX : ufs.particles(Cuts::pid==_pid || Cuts::pid==20443 || Cuts::pid==445)) {
      	if (XX.children().empty()) continue;
        bool matched = false;
       	map<long,int> nRes = nCount;
       	int ncount = ntotal;
       	findChildren(XX,nRes,ncount);
        for (const Particle & omega : ufs.particles(Cuts::pid==223)) {
          if (omega.parents()[0].pid()==_pid || omega.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(omega,nRes2,ncount2);
          matched = true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if     (XX.pid()==20443) _h[1]->fill(_ecms);
            else if(XX.pid()==  445) _h[2]->fill(_ecms);
            // check decay mode
            else {
              int nOut=0;
              map<long,int> nDecay;
              findChildren2(XX, nDecay,nOut);
              if(nOut==3 && nDecay[443]==1 && nDecay[211]==1 && nDecay[-211]==1)
                _h[0]->fill(_ecms);
            }
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    BinnedHistoPtr<string> _h[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2779452);

}
