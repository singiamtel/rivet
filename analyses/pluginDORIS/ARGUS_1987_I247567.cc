// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief cross section for e+e- -> 2pi+ 2pi- pi0 and rho0 omega
  class ARGUS_1987_I247567 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1987_I247567);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==113 || Cuts::pid==223), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_c[ix], 1+ix, 1, 1);
        for (const string& en : _c[ix].binning().edges<0>()) {
          const size_t idx = en.find("-");
          if(idx!=std::string::npos) {
            const double emin = std::stod(en.substr(0,idx));
            const double emax = std::stod(en.substr(idx+1,string::npos));
            if(inRange(sqrtS()/GeV, emin, emax)) {
              _ecms[ix] = en;
              break;
            }
          }
          else {
            const double end = std::stod(en)*GeV;
            if (isCompatibleWithSqrtS(end)) {
              _ecms[ix] = en;
              break;
            }
          }
        }
      }
      if(_ecms[0].empty() && _ecms[1].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int & ncount) {
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      if (ntotal==5 && nCount[211]==2 && nCount[-211]==2 && nCount[111]==1) _c[0]->fill(_ecms[0]);
      // rho0 mesons
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& rho : ufs.particles(Cuts::pid==113)) {
        if (rho.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(rho,nRes,ncount);
        bool matched=false;
        // omega
        for (const Particle& omega : ufs.particles(Cuts::pid==223)) {
       	  if(omega.children().empty()) continue;
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(omega,nRes2,ncount2);
          if(ncount2 !=0 ) continue;
          matched = true;
      	  for (const auto& val : nRes2) {
      	    if (val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
          if (matched) {
            _c[1]->fill(_ecms[1]);
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c, crossSection()/nanobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c[2];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1987_I247567);

}
