// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for D+D- pi+ pi- (D1 D and psi(3770) pi+pi- subprocesses)
  class BESIII_2019_I1756876 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1756876);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_nD1D, 1, 1,  1);
      book(_nPsi, 1, 1,  2);
      for (const string& en : _nD1D.binning().edges<0>()) {
        const double end = std::stod(en)*MeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
        if(child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else
          findChildren(child,nRes,ncount);
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
      for(const Particle& p1 : ufs.particles(Cuts::abspid==411)) {
        int sign = p1.pid()/411;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        matched=false;
        for(const Particle& p2 : ufs.particles(Cuts::pid==-sign*411)) {
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if(ncount2!=2) continue;
          matched=true;
          for(auto const & val : nRes2) {
            if(abs(val.first)==211) {
              if(val.second!=1) {
                matched = false;
                break;
              }
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) break;
        }
        if(matched)
          break;
      }
      if(!matched) vetoEvent;
      // psi(3770) pi+pi-
      for(const Particle& p1 : ufs.particles(Cuts::abspid==30443)) {
        if(p1.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        // psi(3770) pi+pi-
        if(ncount!=2) continue;
        bool matched = true;
        for(auto const & val : nRes) {
          if(abs(val.first)==211) {
            if(val.second !=1) {
              matched = false;
              break;
            }
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if(matched) {
          _nPsi->fill(_ecms);
          return;
        }
      }
      // D1 D
      for(const Particle& p1 : ufs.particles(Cuts::abspid==10413)) {
        if(p1.children().empty()) continue;
        int sign = p1.pid()/10413;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        matched=false;
        for(const Particle& p2 : ufs.particles(Cuts::pid==-sign*411)) {
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if(ncount2!=0) continue;
          matched=true;
          for(auto const & val : nRes2) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) break;
        }
        if(matched) {
          _nD1D->fill(_ecms);
          return;
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/picobarn;
      scale(_nD1D,fact);
      scale(_nPsi,fact);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nD1D, _nPsi;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1756876);

}
