// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross-section for e+e- -> pi0pi0psi(2S) for energies between 4.009 to 4.600 GeV
  class BESIII_2018_I1633425 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1633425);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nPsi2, 1, 1, 1);
      for (const string& en : _nPsi2.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        // find the psi(2S)
        if(p.pid()==100443) {
          map<long,int> nRes = nCount;
          int ncount = ntotal;
          findChildren(p,nRes,ncount);
          // omega pi+pi-
          if(ncount!=2) continue;
          bool matched = true;
          for(auto const & val : nRes) {
            if(abs(val.first)==111) {
              if(val.second !=2) {
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
            _nPsi2->fill(_ecms);
            break;
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nPsi2, crossSection()/ sumOfWeights() /picobarn);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nPsi2;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1633425);

}
