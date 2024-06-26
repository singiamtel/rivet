// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- -> omega chi_b(0,1,2)
  class BELLE_2022_I2142648 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2142648);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
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
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int &ncount) {
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
      for (const Particle & chi : ufs.particles(Cuts::pid==10551 || Cuts::pid==20553 ||  Cuts::pid==555)) {
        bool matched = false;
        if (chi.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(chi,nRes,ncount);
        // loop over omega mesons
        for (const Particle& omega : ufs.particles(Cuts::pid==223)) {
          Particle parent = omega;
          while (!parent.parents().empty()) {
            parent = parent.parents()[0];
            if (parent.pid()==555 || parent.pid()==20553 || parent.pid()==10551) break;
          }
          if ((parent.pid()==555 || parent.pid()==20553 || parent.pid()==10551) &&
              fuzzyEquals(parent.mom(),chi.mom())) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(omega,nRes2,ncount2);
          matched = true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (!matched) continue;
          if (chi.pid()==10551)     _sigma[0]->fill(_ecms);
          else if (chi.pid()==20553) _sigma[1]->fill(_ecms);
          else if (chi.pid()==555)   _sigma[2]->fill(_ecms);
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
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2142648);

}
