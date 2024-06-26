// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> p KS0 nbar K-
  class BESIII_2018_I1681638 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1681638);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==102134), "UFS");
      // counters
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_sigma[ix], 1, 1, 1+ix);
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
      // find the final-state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // FS
      if (ntotal==4&&nCount[310]==1) {
        if ((nCount[ 2212]==1 && nCount[-2112]==1 && nCount[-321]==1) ||
           (nCount[-2212]==1 && nCount[ 2112]==1 && nCount[ 321]==1)) {
          _sigma[0]->fill(_ecms);
        }
      }
      for (const Particle& p : apply<FinalState>(event, "UFS").particles()) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if (ncount!=2) continue;
        bool matched = true;
        int sign = p.pid()/p.abspid();
        for (const auto& val : nRes) {
          if (val.first==-sign*2112 || val.first==310) {
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
          _sigma[1]->fill(_ecms);
          break;
        }
        matched = true;
        for (const auto& val : nRes) {
          if(val.first==-sign*2212 ||
             val.first==sign*321) {
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
          _sigma[2]->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/picobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1681638);

}
