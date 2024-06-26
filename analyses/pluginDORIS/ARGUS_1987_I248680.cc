// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> K*0K*0
  class ARGUS_1987_I248680 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1987_I248680);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV,1.6,3.5)) {
        for (unsigned int ix=0; ix<9; ++ix) {
          book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
        }
      }
      else {
        throw Error("Invalid CMS energy for ARGUS_1987_I248680");
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle& child : p.children()) {
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find any K* mesons
      int ires=-1;
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles Kstar=ufs.particles(Cuts::abspid==313);
      for (unsigned int ix=0;ix<Kstar.size();++ix) {
       	if(Kstar[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(Kstar[ix],nRes,ncount);
        bool matched=false;
        // K*K*
        for (unsigned int iy=ix+1; iy<Kstar.size(); ++iy) {
          if (Kstar[iy].children().empty()) continue;
          if (Kstar[ix].pid()!=-Kstar[iy].pid()) continue;
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(Kstar[iy],nRes2,ncount2);
          if (ncount2 !=0 ) continue;
          matched = true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            break;
          }
        }
        if (matched) {
          _nMeson[1]->fill();
          ires=7;
          break;
        }
        int sign = Kstar[ix].pid()/Kstar[ix].abspid();
        // three body intermediate states
        if (ncount==2) {
          // K*0 K- pi+ +ccd
          matched=true;
          for (const auto& val : nRes) {
            if (val.first==sign*211 || val.first==-sign*321) {
              if (val.second!=1) {
                matched = false;
                break;
              }
            }
            else {
              if (val.second!=0) {
                matched = false;
                break;
              }
            }
          }
          if(matched) {
            _nMeson[2]->fill();
            ires=6;
            break;
          }
        }
      }
      // look for phi modes
      for (const Particle& p : ufs.particles(Cuts::pid==PID::PHI)) {
       	if (p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
        if (ncount==2) {
          bool matched=true;
          for (const auto& val : nRes) {
            if (abs(val.first)==211) {
              if (val.second!=1) {
                matched = false;
                break;
              }
            }
            else {
              if (val.second!=0) {
                matched = false;
                break;
              }
            }
          }
          if (matched) {
            ires=8;
            break;
          }
        }
      }
      // 4 meson modes
      if (ntotal==4 &&
         nCount[PID::KPLUS ]==1 && nCount[PID::KMINUS ]==1 &&
         nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 ) {
          _nMeson[0]->fill();
          _nMeson[4]->fill();
        if (ires<0) {
          _nMeson[3]->fill();
          _nMeson[5]->fill();
        }
        else _nMeson[ires]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nMeson, crossSection()/nanobarn/sumOfWeights());
      // loop over tables in paper
      for (unsigned int ih=1; ih<5; ++ih) {
       	for (unsigned int iy=1; iy<2; ++iy) {
          unsigned int iloc=ih+iy-2;
       	  assert(iloc<=8);
       	  BinnedEstimatePtr<string> mult;
       	  book(mult, ih, 1, iy);
          for (auto& b : mult->bins()) {
            if (isCompatibleWithSqrtS(std::stod(b.xEdge())/GeV)) {
              b.setVal(_nMeson[iloc]->val());
              b.setErr(_nMeson[iloc]->err());
            }
          }
        }
      }
      for (unsigned int iy=1; iy<6; ++iy) {
        unsigned int iloc=3+iy;
        assert(iloc<=8);
        Estimate1DPtr mult;
        book(mult, 5, 1, iy);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS(), b.xMin(), b.xMax())) {
            b.setVal(_nMeson[iloc]->val());
            b.setErr(_nMeson[iloc]->err());
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1987_I248680);

}
