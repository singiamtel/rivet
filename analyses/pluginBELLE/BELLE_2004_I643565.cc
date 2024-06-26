// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > D(*)+ D(*)-
  class BELLE_2004_I643565 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I643565);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==411 || Cuts::abspid==413), "UFS");
      // histos
      book(_h_sigma[0],1,1,2);
      book(_h_sigma[1],1,1,5);
      book(_h_hel     ,2,1,1);
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      Particle Dstar;
      for (unsigned int ix=0;ix<ufs.particles().size();++ix) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(ufs.particles()[ix],nRes,ncount);
        bool matched=false;
        for (unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
          if (ufs.particles()[ix].pid()*ufs.particles()[iy].pid()>0) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(ufs.particles()[iy],nRes2,ncount2);
          if (ncount2!=0) continue;
          matched=true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            if (ufs.particles()[ix].abspid()==413 &&
                ufs.particles()[iy].abspid()==413) {
              _h_sigma[0]->fill("10.58"s);
            }
            else if (ufs.particles()[ix].abspid()==411 &&
                     ufs.particles()[iy].abspid()==413) {
              _h_sigma[1]->fill("10.58"s);
              Dstar = ufs.particles()[iy];
            }
            else if (ufs.particles()[iy].abspid()==411 &&
                     ufs.particles()[ix].abspid()==413) {
              _h_sigma[1]->fill("10.58"s);
              Dstar = ufs.particles()[ix];
            }
            break;
          }
        }
        if (matched) break;
      }
      if (Dstar.abspid()==413 && Dstar.children().size()==2) {
        Particle pim;
        int sign = Dstar.pid()/Dstar.abspid();
        if (Dstar.children()[0].pid()==sign*211 &&
            Dstar.children()[1].pid()==sign*421) {
          pim = Dstar.children()[0];
        }
        else if (Dstar.children()[1].pid()==sign*211 &&
          Dstar.children()[0].pid()==sign*421) {
          pim = Dstar.children()[1];
        }
        else {
          return;
        }
        const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(Dstar.mom().betaVec());
        FourMomentum pPim = boost2.transform(pim.mom());
        _h_hel->fill(pPim.p3().unit().dot(Dstar.mom().p3().unit()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma, crossSection()/picobarn/sumOfWeights());
      normalize(_h_hel, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_sigma[2];
    Histo1DPtr _h_hel;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2004_I643565);

}
