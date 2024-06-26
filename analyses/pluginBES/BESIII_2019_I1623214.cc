// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > eta Y(2175)
  class BESIII_2019_I1623214 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1623214);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 30333);
      // projections
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
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
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
      // loop over any eta mesons
      for (const Particle & eta : ufs.particles(Cuts::pid==221)) {
        bool matched = false;
        if(eta.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(eta,nRes,ncount);
        for (const Particle & phi : ufs.particles(Cuts::pid==333)) {
          if (phi.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(phi,nRes2,ncount2);
          matched = true;
          // required eta phi pi+ pi- final state
          for (auto const & val : nRes2) {
            if (abs(val.first)==211) {
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
          if (!matched) continue;
          matched = false;
          // finally check we have f0
          for (const Particle & f0 : ufs.particles(Cuts::pid==9010221)) {
            if (f0.children().empty()) continue;
            map<long,int> nRes3 = nRes2;
            int ncount3 = ncount2;
            findChildren(f0,nRes3,ncount3);
            matched = true;
            for (const auto& val : nRes3) {
              if (val.second!=0) {
                matched = false;
                break;
              }
            }
            if (matched) break;
          }
          if (matched) break;
        }
        if (!matched) continue;
        // finally check phi(3D1) present
        for (const Particle & phi : ufs.particles(Cuts::pid==_pid)) {
          if (phi.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(phi,nRes2,ncount2);
          matched = true;
          // required eta phi(3D1) final state
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) break;
        }
        if (matched) {
          _sigma->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma,crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    BinnedHistoPtr<string> _sigma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1623214);

}
