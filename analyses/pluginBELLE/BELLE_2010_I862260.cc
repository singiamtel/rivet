// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> eta eta
  class BELLE_2010_I862260 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2010_I862260);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV,1.096,4.)) {
        book(_sigmaEtaEta[0],"TMP/nEtaPi_1",refData(1, 1, 1));
        if (sqrtS()<=2.) book(_sigmaEtaEta[1],"TMP/nEtaPi_2",refData(1, 1, 2));
       	double sMin=1.096,sMax=1.12, step=0.04;
       	unsigned int ihist=2;
       	while (sMin<3.3) {
          if (inRange(sqrtS()/GeV, sMin, sMax)) {
            break;
          }
          sMin=sMax;
          sMax+=step;
          ihist+=1;
          if (fuzzyEquals(2.4, sMin)) step=0.1;
        }
        if (ihist<=43) book(_h_cTheta,ihist,1,1);
      }
      else {
        throw Error("Invalid CMS energy for BELLE_2010_I862260");
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
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
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles etas=ufs.particles(Cuts::pid==PID::ETA);
      if (etas.size()<2) vetoEvent;
      const FinalState& fs = apply<FinalState>(event, "FS");
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      for (unsigned int ix=0;ix<etas.size();++ix) {
       	if (etas[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(etas[ix],nRes,ncount);
        for (unsigned int iy=ix+1; iy<etas.size(); ++iy) {
          if (etas[iy].children().empty()) continue;
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(etas[iy],nRes2,ncount2);
          if (ncount2 !=0 ) continue;
          bool matched = true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            double cTheta = abs(etas[iy].momentum().z()/etas[iy].momentum().p3().mod());
            if (cTheta<=0.9) _sigmaEtaEta[0]->fill(sqrtS());
            if (_sigmaEtaEta[1]) _sigmaEtaEta[1]->fill(sqrtS());
            if (_h_cTheta ) _h_cTheta ->fill(cTheta);
            break;
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) scale(_h_cTheta ,fact);
      for (unsigned int ix=0;ix<2;++ix) {
        if (!_sigmaEtaEta[ix]) continue;
        scale(_sigmaEtaEta[ix], fact);
        Estimate1DPtr tmp;
        book(tmp, 1, 1, 1+ix);
        barchart(_sigmaEtaEta[ix],tmp);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaEtaEta[2], _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2010_I862260);

}
