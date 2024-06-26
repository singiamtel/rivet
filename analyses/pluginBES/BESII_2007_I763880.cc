// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross-sections for light hadrons at 3.773 and 3.650 GeV
  class BESII_2007_I763880 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2007_I763880);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=4;ix<8;++ix) {
        book(_nMeson[ix], 1,1,ix);
      }

      for (const string& en : _nMeson[5].binning().edges<0>()) {
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
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        if(p.pid()!=223 && p.pid()!=113&&abs(p.pid())!=313&& abs(p.pid())!=323)
          continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if(p.pid()==113) {
          if(ncount!=3) continue;
          unsigned int npi(0),nK(0);
          bool matched = true;
          for(auto const & val : nRes) {
            if(abs(val.first)==111 && val.second==1) {
              npi+=1;
            }
            else if(abs(val.first)==321 && val.second==1) {
              nK+=1;
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if(npi==1&&nK==2)
              if(_ecms=="3.65") _nMeson[4]->fill(_ecms);
          }
        }
        else if(abs(p.pid())==213) {
          if(ncount!=3) continue;
          unsigned int npi(0),nK(0);
          bool matched = true;
          int ipi = p.pid()==213 ? -211 : 211;
          for(auto const & val : nRes) {
            if(abs(val.first)== ipi && val.second==1) {
              npi+=1;
            }
            else if(abs(val.first)==321 && val.second==1) {
              nK+=1;
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if(npi==1&&nK==2)
              _nMeson[5]->fill(_ecms);
          }
        }
        else if(abs(p.pid())==313) {
          if(ncount!=3) continue;
          unsigned int npi(0),nK(0),npi0(0);
          bool matched = true;
          int ipi = p.pid()==313 ?  211 : -211;
          int iK  = p.pid()==313 ? -321 :  321;
          for(auto const & val : nRes) {
            if(abs(val.first)== ipi && val.second==1) {
              npi+=1;
            }
            else if(abs(val.first)==iK && val.second==1) {
              nK+=1;
            }
            if(abs(val.first)== 111 && val.second==1) {
              npi0+=1;
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if(npi==1&&nK==1&&npi0==1)
              _nMeson[6]->fill(_ecms);
          }
        }
        else if(abs(p.pid())==323) {
          if(ncount!=3) continue;
          unsigned int npi(0),nK(0),npi0(0);
          bool matched = true;
          int ipi = p.pid()==323 ?  211 : -211;
          int iK  = p.pid()==323 ? -321 :  321;
          for(auto const & val : nRes) {
            if(abs(val.first)== ipi && val.second==1) {
              npi+=1;
            }
            else if(abs(val.first)==iK && val.second==1) {
              nK+=1;
            }
            if(abs(val.first)== 111 && val.second==1) {
              npi0+=1;
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if(npi==1&&nK==1&&npi0==1)
              _nMeson[7]->fill(_ecms);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn; 
      for(unsigned int ix=4;ix<8;++ix)
        scale(_nMeson[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nMeson[8];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2007_I763880);

}
