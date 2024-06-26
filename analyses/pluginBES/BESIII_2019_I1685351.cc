// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+ e- > mu+ mu-, pi+ pi- eta and 5pi near J/psi
  class BESIII_2019_I1685351 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1685351);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==PID::ETA), "UFS");
      // histograms
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_sigma[ix], 1, 1, 1+ix);
      }
      // central beam energy
      string ecms = std::to_string(0.01*double(round(sqrtS()/MeV*100)));
      const size_t idx = ecms.find(".");
      ecms = ecms.substr(0,idx+3);
      if(ecms[ecms.length()-1]=='0') ecms.pop_back();
      _eCent = getOption<string>("ECENT", ecms);
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
      // mu+mu- + photons
      if (nCount[-13]==1 && nCount[13]==1 && ntotal==2+nCount[22]) {
        _sigma[0]->fill(_eCent);
      }
      else if (nCount[111]==1 && nCount[211]==2 && nCount[-211]==2  && ntotal==5+nCount[22]) {
        _sigma[2]->fill(_eCent);
      }

      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over eta mesons
      for (const Particle& p : ufs.particles()) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==211) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if (val.first!=PID::PHOTON && val.second!=0) {
            matched = false;
            break;
          }
        }
        if (!matched) continue;
        _sigma[1]->fill(_eCent);
        break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _eCent;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1685351);

}
