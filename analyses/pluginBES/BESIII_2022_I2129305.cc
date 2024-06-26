// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > pi+pi- D+ D-
  class BESIII_2022_I2129305 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2129305);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0;ix<2;++ix) {
        book(_sigma[ix], ix+1, 1, 1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
    }

    void findChildren(const Particle& p,map<long,int> & nRes, int &ncount) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool matched=false;
      for (const Particle& Dp : ufs.particles(Cuts::pid==PID::DPLUS)) {
        if (Dp.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(Dp,nRes,ncount);
        for (const Particle& Dm : ufs.particles(Cuts::pid==PID::DMINUS)) {
          if (Dm.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(Dm,nRes2,ncount2);
          if (ncount2!=2) continue;
          matched=true;
          for (const auto& val : nRes2) {
            if (abs(val.first)==211) {
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
          if(matched) {
            _sigma[0]->fill(_ecms[0]);
            if (!Dm.parents().empty() && !Dp.parents().empty()) {
              // see if D+ from psi3
              Particle parent1=Dp;
              while (!parent1.parents().empty()) {
                parent1=parent1.parents()[0];
                if (parent1.pid()==447) break;
              }
              if (parent1.pid()!=447) break;
              // see if D- from psi3
              Particle parent2=Dm;
              while (!parent2.parents().empty()) {
                parent2=parent2.parents()[0];
                if (parent2.pid()==447) break;
              }
              if (parent2.pid()!=447) break;
              if (fuzzyEquals(parent1.momentum(),parent2.momentum())) {
                _sigma[1]->fill(_ecms[1]);
              }
            }
            break;
          }
        }
        if(matched) break;
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
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2129305);

}
