// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- -> gamma X(3872)
  class BESIII_2019_I1724671 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1724671);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counters
      for (unsigned int ix=0;ix<2;++ix) {
        book(_sigma[ix],1,1,1+ix);
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

    void findChildren(const Particle & p,map<long,int>& nRes, int &ncount) {
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
      int mode=-1;
      for (const Particle & psi : ufs.particles(Cuts::pid==443)) {
        if (psi.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(psi,nRes,ncount);
        // check for J/psi pi+pi- mode
        if (ncount==3) {
          bool matched = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==211 || val.first==PID::GAMMA) {
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
            mode=1;
            break;
          }
        }
        bool matched = false;
        for (const Particle& omega : ufs.particles(Cuts::pid==223)) {
          if (omega.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(omega,nRes2,ncount2);
          if (ncount2!=1) continue;
          matched = true;
          for (const auto& val : nRes2) {
            if (val.first==PID::GAMMA) {
              if (val.second !=1) {
                matched = false;
                break;
              }
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) break;
        }
        if (matched) {
          mode = 0;
          break;
        }
      }
      if (mode<0) vetoEvent;
      // find the intermediate
      for (const Particle& XX : ufs.particles(Cuts::pid==_pid)) {
        if (XX.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(XX, nRes, ncount);
        if (ncount!=1) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (val.first==PID::GAMMA) {
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
          _sigma[mode]->fill(_ecms);
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
    int _pid;
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1724671);

}
