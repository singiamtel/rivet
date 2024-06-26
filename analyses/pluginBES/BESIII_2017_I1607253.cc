// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > 2 phi + phi or omega
  class BESIII_2017_I1607253 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1607253);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==223 || Cuts::pid==333), "UFS");
      // counters
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_sigma[ix], 1+ix, 1, 1);
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
      for (const Particle& p: fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      Particles phi   = ufs.particles(Cuts::pid==333);
      Particles omega = ufs.particles(Cuts::pid==223);
      if (phi.size()<2) vetoEvent;
      bool found = false;
      // first phi
      for (unsigned int ix=0; ix<phi.size(); ++ix) {
        const Particle & phi1 = phi[ix];
        if(phi1.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(phi1,nRes,ncount);
        // second phi
        for (unsigned int iy=ix+1; iy<phi.size(); ++iy) {
          const Particle & phi2 = phi[iy];
          if (phi2.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(phi2,nRes2,ncount2);
          // search for third phi
          for (unsigned int iz=iy+1; iz<phi.size(); ++iz) {
            const Particle & phi3 = phi[iz];
            if (phi3.children().empty()) continue;
            map<long,int> nRes3 = nRes2;
            int ncount3 = ncount2;
            findChildren(phi3,nRes3,ncount3);
            found = true;
            for (const auto& val : nRes3) {
              if (val.second!=0) {
                found = false;
                break;
              }
            }
            if(found) {
              _sigma[1]->fill(_ecms);
              break;
            }
          }
          if (found) break;
          // search for omega
          for (unsigned int iz=0; iz<omega.size(); ++iz) {
            const Particle& omega1 = omega[iz];
            if (omega1.children().empty()) continue;
            map<long,int> nRes3 = nRes2;
            int ncount3 = ncount2;
            findChildren(omega1,nRes3,ncount3);
            found = true;
            for (const auto& val : nRes3) {
              if (val.second!=0) {
                found = false;
                break;
              }
            }
            if (found) {
              _sigma[0]->fill(_ecms);
              break;
            }
          }
          if(found) break;
        }
        if (found) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /femtobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1607253);

}
