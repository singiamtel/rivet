// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > pi0 pi0 hc
  class BESIII_2014_I1318650 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2014_I1318650);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_sigma[ix], 1, 1, 1+ix);
      }

      for (const string& en : _sigma[0].binning().edges<0>()) {
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
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      // loop over any h_c
      for (const Particle& hc : ufs.particles(Cuts::pid==10443)) {
        if (hc.children().empty()) continue;
        // find the h_c
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(hc,nRes,ncount);
        // h_c pi0 pi0
        if(ncount!=2) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==111) {
            if (val.second !=2) {
              matched = false;
              break;
            }
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
          _sigma[0]->fill(_ecms);
          if (hc.parents().empty()) break;
          Particle Zc = hc.parents()[0];
          if (Zc.pid()==_pid && Zc.children().size()==2 &&
             (Zc.children()[0].pid()==PID::PI0 ||
              Zc.children()[1].pid()==PID::PI0)) _sigma[1]->fill(_ecms);
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
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2014_I1318650);

}
