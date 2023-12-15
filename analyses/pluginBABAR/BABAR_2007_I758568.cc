// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2007_I758568 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I758568);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=1;ix<12;++ix) {
        stringstream ss;
        ss << "TMP/n" << ix;
        book(_nMeson[ix], ss.str());
      }
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
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        if(p.pid()!=223 && p.pid()!=221&&p.pid()!=331&& p.pid()!=20223)
          continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // omega
        if(p.pid()==223) {
          if(ncount==2) {
            bool matched = true;
            for(auto const & val : nRes) {
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
            if(matched) _nMeson[3]->fill();

            matched = true;
            for(auto const & val : nRes) {
              if(abs(val.first)==321) {
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
            if(matched) _nMeson[10]->fill();
          }

          for (const Particle& p2 : ufs.particles()) {
            if(p2.pid()!=9010221) continue;
            if(p2.parents()[0].isSame(p)) continue;
            map<long,int> nResB = nRes;
            int ncountB = ncount;
            findChildren(p2,nResB,ncountB);
            if(ncountB!=0) continue;
            bool matched2 = true;
            for(auto const & val : nResB) {
              if(val.second!=0) {
                matched2 = false;
                break;
              }
            }
            if(matched2) _nMeson[4]->fill();
          }
        }
        else if(p.pid()==221) {
          if(ncount==2) {
            bool matched = true;
            for(auto const & val : nRes) {
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
            if(matched) _nMeson[2]->fill();
          }
          else if(ncount==4) {
            bool matched = true;
            for(auto const & val : nRes) {
              if(abs(val.first)==211) {
                if(val.second!=2) {
                  matched = false;
                  break;
                }
              }
              else if(val.second!=0) {
                matched = false;
                break;
              }
            }
            if(matched) _nMeson[5]->fill();

            matched=true;
            for(auto const & val : nRes) {
              if(abs(val.first)==211 || abs(val.first)==321) {
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
            if(matched) _nMeson[11]->fill();
          }

          for(const Particle& p2 : ufs.particles()) {
            if(p2.pid()!=333) continue;
            map<long,int> nResB = nRes;
            int ncountB = ncount;
            findChildren(p2,nResB,ncountB);
            if(ncountB!=0) continue;
            bool matched2 = true;
            for(auto const & val : nResB) {
              if(val.second!=0) {
                matched2 = false;
                break;
              }
            }
            if(matched2) _nMeson[9]->fill();
          }
        }
        else if(p.pid()==331 || p.pid()==20223) {
          if(ncount!=2) continue;
          bool matched = true;
          for(auto const & val : nRes) {
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
          if(matched) {
            if(p.pid()==331)
              _nMeson[6]->fill();
            else
              _nMeson[7]->fill();
          }
        }
      }

      if(ntotal==5) {
        if(nCount[211]==2&&nCount[-211]==2&&nCount[111]==1)
          _nMeson[1]->fill();
        else if(nCount[321]==1&&nCount[-321]==1&&
                nCount[211]==1&&nCount[-211]==1&&nCount[111]==1)
          _nMeson[8]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<12;++ix) {
        double sigma = _nMeson[ix]->val();
        double error = _nMeson[ix]->err();
        sigma *= crossSection()/ sumOfWeights() /nanobarn;
        error *= crossSection()/ sumOfWeights() /nanobarn;
        Estimate1DPtr  mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[12];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I758568);


}
