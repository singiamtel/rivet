// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > pi+ pi- h_c
  class BESIII_2017_I1494065 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1494065);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==10443), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
        book(_nhc[ix], 1+ix, 1, 1);
        for (const string& en : _nhc[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over the h_c
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // h_c pi+pi-
        if(ncount!=2) continue;
        bool matched = true;
        for(auto const & val : nRes) {
          if(abs(val.first)==211) {
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
        if(matched) {
          _nhc[0]->fill(_ecms[0]);
          _nhc[1]->fill(_ecms[1]);
          break;
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<2;++ix)
        scale(_nhc[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nhc[2];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1494065);

}
