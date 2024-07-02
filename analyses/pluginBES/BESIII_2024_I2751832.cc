// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- -> pi+ pi- pi0
  class BESIII_2024_I2751832 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2751832);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix)
	book(_sigma[ix], 1, 1, 1+ix);
      
      for (const string& en : _sigma[0].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int& ncount) {
      for (const Particle& child : p.children()) {
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
      Particle pip, pim;
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
        if (p.pid()     == 211) pip=p;
        else if (p.pid()==-211) pim=p;
      }
      if (ntotal!=3) vetoEvent;
      if (nCount[-211]==1&&nCount[211]==1&&nCount[111]==1) {
        _sigma[0]->fill(_ecms);
      }
      else {
        vetoEvent;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==113 || Cuts::abspid==213 ||
                                             Cuts::abspid==100113 || Cuts::abspid==100213)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if (ncount!=1) continue;
        int idOther = 211;
        if (p.pid()==113 || p.pid()==100113) {
          idOther = 111;
        }
        else if (p.pid()==213 || p.pid()==100213) {
          idOther = -211;
        }
        bool matched=true;
        for (const auto& val : nRes) {
          if (val.first==idOther) {
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
        if (!matched) continue;
        if (matched) {
          if (p.pid()==213 || p.pid()==113)
	    _sigma[1]->fill(_ecms);
	  else
	    _sigma[2]->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix)
	scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2751832);

}
