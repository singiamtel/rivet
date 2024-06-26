// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > omega pi0, K* Kbar K2 Kbar
  class BELLE_2013_I1252555 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1252555);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counters
      for (unsigned int ix=0; ix<5; ++ix) {
        book(_sigma[ix],1+ix,1,1);
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
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // resonances
      for (const Particle & p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==223 ||
										Cuts::abspid==313 || Cuts::abspid==323 ||	Cuts::abspid==315 || Cuts::abspid==325)) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if (ncount!=1) continue;
        unsigned int imode=0;
        long pid = 111;
        if (p.abspid()==313) {
          imode=1;
          pid=310;
        }
        else if (p.abspid()==323) {
          imode=2;
          pid = p.pid()>0 ? -321 : 321;
        }
        else if (p.abspid()==315) {
          imode=3;
          pid=310;
        }
        else if (p.abspid()==325) {
          imode=4;
          pid = p.pid()>0 ? -321 : 321;
        }
        bool matched=true;
        for (const auto& val : nRes) {
          if ((pid==310 && (val.first==310 || val.second==130)) || val.first==pid) {
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
        if (matched) {
          _sigma[imode]->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /femtobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[5];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1252555);

}
