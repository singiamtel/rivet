// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > Delta++ pbar- pi-
  class BESIII_2023_I2661512 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2661512);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==2224), "UFS");
      vector<string> energies({"2.3094", "2.3864", "2.396", "2.5", "2.6444", "2.6464"});
      bool matched=false;
      for(const string& en : energies) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          matched=true;
          break;
        }
      }
      if(!_ecms.empty())
        book(_sigma[0],1,1,1);
      if(isCompatibleWithSqrtS(2.6454)) {
        book(_sigma[1],2,1,1);
        matched=true;
      }
      if(!matched) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int& ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child, nRes, ncount);
        }
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
      // loop over delta++
      for (const Particle& p : apply<FinalState>(event, "UFS").particles()) {
        int sign =-p.pid()/p.abspid();
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p, nRes, ncount);
        bool matched = true;
        for (const auto& val : nRes) {
          if (val.first== 211*sign || val.first==2212*sign) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if (val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
	  if(_sigma[0]) _sigma[0]->fill(_ecms);
          if(_sigma[1]) _sigma[1]->fill(string("2.6454"));
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<2;++ix) scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2661512);

}
