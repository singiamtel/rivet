// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> KK 3pi
  class BABAR_2022_I2120528 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2022_I2120528);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<10; ++ix) {
        book(_sigma[ix], "TMP/n" + toString(ix+1), refData(1+ix,1,1));
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
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
      // stable particles
      if (ntotal==5) {
        // K+K- 3pi0
        if(nCount[-321]==1 && nCount[321]==1 && nCount[111]==3) _sigma[0]->fill(sqrtS());
        // KS0 Kpi 2pi0
        else if (nCount[111]==2 &&nCount[310]==1 &&
               ((nCount[ 321]==1 && nCount[-211]==1) ||
                (nCount[-321]==1 && nCount[ 211]==1))) _sigma[1]->fill(sqrtS());
        // KS0 Kpi pi+pi-
        else if (nCount[310]==1 &&
               ((nCount[ 321]==1 && nCount[-211]==2 && nCount[ 211]==1) ||
                (nCount[-321]==1 && nCount[ 211]==2 && nCount[-211]==1))) _sigma[2]->fill(sqrtS());
      }
      // now the resonant states
      // first K+K- eta and phi eta
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==PID::ETA)) {
        if (p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched1  = false;
        if (ncount==2) {
          matched1  = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==PID::KPLUS) {
              if (val.second!=1) {
                matched1 = false;
                break;
              }
            }
            else if (val.second!=0) {
              matched1 = false;
              break;
            }
          }
          if (matched1) {
            _sigma[3]->fill(sqrtS());
          }
        }
        bool matched2=false;
        for (const Particle& p2 : ufs.particles(Cuts::pid==PID::PHI)) {
          if (p2.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2!=0) continue;
          matched2=true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched2 = false;
              break;
            }
          }
          if (matched2) {
            _sigma[4]->fill(sqrtS());
            break;
          }
        }
        if (matched1 || matched2) break;
      }
      // states with a K*
      Particles Kstar = ufs.particles(Cuts::abspid==313 or Cuts::abspid==323);
      for (unsigned int ix=0; ix<Kstar.size(); ++ix) {
        if (Kstar[ix].children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(Kstar[ix],nRes,ncount);
        bool matched1=false;
        if (ncount==3) {
          int idother=310;
          if (Kstar[ix].abspid()==323) idother= -321*(Kstar[ix].pid()/Kstar[ix].abspid());
          for (const auto& val : nRes) {
            if (val.first==idother) {
              if (val.second!=1) {
                matched1 = false;
                break;
              }
            }
            else if (val.first==PID::PI0) {
              if(val.second!=2) {
                matched1 = false;
                break;
              }
            }
            else if(val.second!=0) {
              matched1 = false;
              break;
            }
          }
          if (matched1) {
            if (Kstar[ix].abspid()==321) _sigma[7]->fill(sqrtS());
            else                         _sigma[6]->fill(sqrtS());
          }
          // K*+ K*-
          bool matched2=false;
          for (unsigned int iy=0;iy<ix+1;++iy) {
            if (Kstar[iy].children().empty() ||
                Kstar[iy].abspid()==313 ||
                Kstar[iy].pid()!=-Kstar[ix].pid()) {
              continue;
            }
            map<long,int> nRes2 = nRes;
            int ncount2 = ncount;
            findChildren(Kstar[iy],nRes2,ncount2);
            if (ncount2!=1) continue;
              matched2=true;
              for (const auto& val : nRes2) {
                if (val.first==PID::PI0) {
                  if (val.second!=1) {
                    matched2 = false;
                    break;
                  }
                }
                if (val.second!=0) {
                  matched2 = false;
                  break;
                }
              }
              if (matched2) {
                _sigma[5]->fill(sqrtS());
                break;
              }
            }
          }
        }
        // finally f_1
        for (const Particle& p : ufs.particles(Cuts::abspid==PID::RHOPLUS)) {
          if (p.children().empty()) continue;
          map<long,int> nRes=nCount;
          int ncount = ntotal;
          findChildren(p,nRes,ncount);
          bool matched1  = false;
          if (ncount==3) {
            int idother = -(p.pid()/p.abspid())*321;
            matched1  = true;
          for(const auto& val : nRes) {
            if (val.first==PID::PI0 || val.first==idother || val.first==PID::K0S) {
              if (val.second!=1) {
                matched1 = false;
                break;
              }
            }
            else if (val.second!=0) {
              matched1 = false;
              break;
            }
          }
          if (matched1) {
            _sigma[8]->fill(sqrtS());
          }
        }
        if (matched1) break;
      }
      // finally f_1
      for (const Particle& p : ufs.particles(Cuts::pid==20223)) {
        if (p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched1  = false;
        if (ncount==2) {
          matched1  = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==PID::PIPLUS) {
              if (val.second!=1) {
                matched1 = false;
                break;
              }
            }
            else if (val.second!=0) {
              matched1 = false;
              break;
            }
          }
          if (matched1) {
            _sigma[9]->fill(sqrtS());
          }
        }
        if (matched1) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<10;++ix) {
        scale(_sigma[ix],fact);
        Estimate1DPtr tmp;
        book(tmp,1+ix,1,1);
        barchart(_sigma[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[10];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2022_I2120528);

}
