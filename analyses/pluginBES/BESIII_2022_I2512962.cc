// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> eta lambda lambda bar
  class BESIII_2022_I2512962 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2512962);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counter
      book(_sigma,1,1,1);

      for (const string& en : _sigma.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over eta mesons
      for (const Particle& eta : ufs.particles(Cuts::pid==221)) {
        bool matched = false;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(eta,nRes,ncount);
        // then Lambda baryons
        for (const Particle& lambda : ufs.particles(Cuts::pid==3122)) {
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(lambda,nResB,ncountB);
          for (const Particle& lambar : ufs.particles(Cuts::pid==-3122)) {
            map<long,int> nResC = nResB;
            int ncountC = ncountB;
            findChildren(lambar,nResC,ncountC);
            matched=true;
            for (const auto& val : nResC) {
              if (val.second!=0) {
                matched = false;
                break;
              }
            }
            if (matched) {
              _sigma->fill(_ecms);
              break;
            }
          }
          if (matched) break;
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma,crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2512962);

}
