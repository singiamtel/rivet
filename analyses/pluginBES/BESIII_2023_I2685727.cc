// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > phi K+K-/KS0 KS0
  class BESIII_2023_I2685727 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2685727);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<2; ++ix)
        book(_nPhi[ix], 1, 1, 1+ix);
      for (const string& en : _nPhi[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
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
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // find the phis
      for (const Particle& p : ufs.particles(Cuts::pid==333)) {
      	if (p.children().empty()) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p,nRes,ncount);
      	// phi KK
        if (ncount!=2) continue;
        bool matched = false;
        for (unsigned int imode=0; imode<2; ++imode) {
          matched = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==310 && imode==1) {
              if (val.second !=2) {
                matched = false;
                break;
              }
            }
            else if (abs(val.first)==321 && imode==0) {
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
            _nPhi[imode]->fill(_ecms);
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nPhi, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nPhi[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2685727);

}
