// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief gamma gamma -> Lambda0 Lambdabar0
  class CLEOII_1997_I439745 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1997_I439745);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(Cuts::abspid==3122 || Cuts::abspid==3212), "UFS");
      // gamma gamma collision
      if (beamIDs().first==22 && beamIDs().second==22)  book(_c,"TMP/nLam",refData(1,1,1));
      else {
        for (unsigned int ix=0; ix<2; ++ix) {
          book(_h_sigma[ix], 2+ix, 1, 1);
        }
      }
    }

    void findChildren(const Particle& p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        } else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    bool findScattered(const Particle& beam, double& q2) {
      bool found = false;
      Particle scat = beam;
      while (!scat.children().empty()) {
        found = false;
        for (const Particle & p : scat.children()) {
          if (p.pid()==scat.pid()) {
            scat=p;
            found=true;
            break;
          }
        }
        if (!found) break;
      }
      if (!found) return false;
      q2 = -(beam.mom() - scat.mom()).mass2();
      return true;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // e+e- collision
      int ne=0;
      const Beam& beams = apply<Beam>(event, "Beams");
      if (_h_sigma[0]) {
        // find scattered leptons and calc Q2
        double q12 = -1, q22 = -1;
        if (!findScattered(beams.beams().first,  q12)) vetoEvent;
        if (!findScattered(beams.beams().second, q22)) vetoEvent;
        ne=1;
      }
      // check the final state
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      Particle Lambda,Lambar;
      // loop over baryon
      bool matched = false;
      for (const Particle& p : ufs.particles()) {
        if (p.children().empty() || p.pid()<0) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // and antibaryons
        for (const Particle& p2 : ufs.particles()) {
          if (p2.children().empty() || p2.pid()>0) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if(ncount2!=2*ne) continue;
          matched=2;
          for (const auto& val : nRes2) {
            if(abs(val.first)==11) {
            if(val.second!=ne) {
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
            Lambda=p;
            Lambar=p2;
            break;
          }
        }
      }
      if (!matched) vetoEvent;
      if (_h_sigma[0]) {
        FourMomentum pCMS = Lambda.mom()+Lambar.mom();
        Vector3 axis = pCMS.p3().unit();
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(pCMS.betaVec());
        if (abs(boost.transform(Lambda.mom()).p3().unit().dot(axis))>0.6) vetoEvent;
        if (abs(boost.transform(Lambar.mom()).p3().unit().dot(axis))>0.6) vetoEvent;
        _h_sigma[0]->fill("10.6"s);
        if (Lambda.pid()==3122 && Lambar.pid()==-3122) _h_sigma[1]->fill("10.6"s);
      }
      else {
        Vector3 axis = beams.beams().first.mom().p3().unit();
        if (abs(axis.dot(Lambda.mom().p3().unit()))>0.6) vetoEvent;
        if (abs(axis.dot(Lambar.mom().p3().unit()))>0.6) vetoEvent;
        if (Lambda.pid()==3122 && Lambar.pid()==-3122) _c->fill(sqrtS()/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_h_sigma[0]) {
        scale(_h_sigma, crossSection()/picobarn/sumOfWeights());
      }
      else {
        scale(_c, crossSection()/nanobarn/sumOfWeights());
        Estimate1DPtr tmp;
        book(tmp,1,1,1);
        barchart(_c,tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _c;
    BinnedHistoPtr<string> _h_sigma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1997_I439745);

}
