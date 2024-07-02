// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- to hadrons
  class BESIII_2015_I1356733 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1356733);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 20445);
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid== 20443 || Cuts::pid==100443), "UFS");
      // histograms
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
      for (const Particle& psi : apply<UnstableParticles>(event, "UFS").particles()) {
        if (psi.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(psi,nRes,ncount);
        int nphoton = psi.pid()==100443 ? 0 : 1;
        if (ncount!=2+nphoton) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==211) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if (val.first==22) {
            if (val.second!=nphoton) {
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
          if (nphoton==0) {
            _sigma[1]->fill(_ecms);
            break;
          }
          else {
            Particle parent = psi.parents()[0];
            if (parent.pid()==_pid && parent.children().size()==2 &&
               (parent.children()[0].pid()==22 || parent.children()[1].pid()==22 )) {
              _sigma[0]->fill(_ecms);
              break;
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
        scale(_sigma[ix], crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1356733);

}
