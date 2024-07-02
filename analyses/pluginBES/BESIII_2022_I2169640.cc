// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- -> phi chi_c(1,2)
  class BESIII_2022_I2169640 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2169640);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_sigma[ix],1,1,1+ix);
      }
      for (const string& en : _sigma[0].binning().edges<0>()) {
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
      // loop over any phi mesons
      for (const Particle & phi : ufs.particles(Cuts::pid==333)) {
        bool matched = false;
        if(phi.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(phi,nRes,ncount);
        for (const Particle & chi : ufs.particles(Cuts::pid==20443 || Cuts::pid==445)) {
          if (chi.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(chi,nRes2,ncount2);
          matched = true;
          for (const auto& val : nRes2) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if (!matched) continue;
          if     (chi.pid()==20443) _sigma[0]->fill(_ecms);
          else if(chi.pid()==  445) _sigma[1]->fill(_ecms);
          break;
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
	    scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2169640);

}
