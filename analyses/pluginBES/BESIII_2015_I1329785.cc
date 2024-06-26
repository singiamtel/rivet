// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to \gamma\chi_{c(0,1,2)}$ at energies between 4.009 and 4.36 GeV
  class BESIII_2015_I1329785 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1329785);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nChi0, 1, 1, 8);
      book(_nChi1, 2, 1, 8);
      book(_nChi2, 3, 1, 8);
      for (const string& en : _nChi0.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
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
      for (const Particle& p :  fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p :  ufs.particles(Cuts::pid==10441 or Cuts::pid==20443 or Cuts::pid==445)) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // chi gamma
        if(ncount!=1) continue;
        bool matched = true;
        for(auto const & val : nRes) {
          if(val.first==22) {
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
          if(p.pid()==10441)
            _nChi0->fill(_ecms);
          else if(p.pid()==20443)
            _nChi1->fill(_ecms);
          else if(p.pid()==445)
            _nChi2->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /picobarn;
      scale(_nChi0,fact);
      scale(_nChi1,fact);
      scale(_nChi2,fact);
    }
    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nChi0,_nChi1,_nChi2;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1329785);


}
