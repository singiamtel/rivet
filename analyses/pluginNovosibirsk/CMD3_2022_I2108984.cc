// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e-> KS0 Kpi pi+pi-
  class CMD3_2022_I2108984 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2022_I2108984);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==20223), "UFS");
      book(_sigma,1,1,1);
      for (const string& en : _sigma.binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/MeV, emin, emax)) {
            _ecms.push_back(en);
          }
        }
        else {
          const double end = std::stod(en)*MeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms.push_back(en);
          }
        }
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_res[ix], "TMP/c_"+toString(ix),refData(2, 1, ix+1));
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      // the FS cross section
      bool foundFS = false;
      if (ntotal==5 && nCount[310]==1 &&
	        ( ( nCount[211]==2 && nCount[-211]==1 && nCount[-321] ==1 ) ||
            ( nCount[211]==1 && nCount[-211]==2 && nCount[ 321] ==1 ) ) ) {
        foundFS = true;
        for(const string & en : _ecms) _sigma->fill(en);
      }
      // check for intermediate f_1
      for (const Particle& p : apply<FinalState>(event, "UFS").particles()) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if(ncount!=2) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==211 ) {
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
          _res[1]->fill(sqrtS()/MeV);
          if (foundFS) _res[0]->fill(sqrtS()/MeV);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      scale(_sigma,fact);
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_res[ix], fact);
        Estimate1DPtr tmp;
        book(tmp,2,1,1+ix);
        barchart(_res[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    Histo1DPtr _res[2];
    vector<string> _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2022_I2108984);

}
