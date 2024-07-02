// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> pi+ pi- pi0 chi_b1,2
  class BELLE_2018_I1678261 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1678261);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==20553 ||	Cuts::pid==555), "UFS");
      book(_sigma,1,1,1);
      for (const string& en : _sigma.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      // loop over any chi mesons
      for (const Particle& chi : ufs.particles()) {
        bool matched = false;
        if (chi.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(chi,nRes,ncount);
        if (ncount==3) {
          matched = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==PID::PIPLUS || val.first==PID::PI0) {
              if (val.second!=1) {
                matched = false;
                break;
              }
            }
            else if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (!matched) continue;
          _sigma->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1678261);

}
