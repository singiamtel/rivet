// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D* rho
  class CLEOII_2003_I611745 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2003_I611745);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 423);
      BB.addStable(-423);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( 213);
      BB.addStable(-213);
      BB.addStable( 113);
      declare(BB, "BB");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<3; ++iy) {
          book(_h[ix][iy], 2, 1+iy, 1+ix);
          if(iy<2) book(_p[ix][iy], "TMP/TMP_"+toString(ix)+"_"+toString(iy));
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
        unsigned int imode=0;
       	if ( ( BB.decaying()[ix].pid() == 521 && BB.modeMatches(ix,2,mode1  ) ) ||
             ( BB.decaying()[ix].pid() ==-521 && BB.modeMatches(ix,2,mode1CC) ) ) imode=0;
       	else if ( ( BB.decaying()[ix].pid() == 511 && BB.modeMatches(ix,2,mode2  ) ) ||
                  ( BB.decaying()[ix].pid() ==-511 && BB.modeMatches(ix,2,mode2CC) ) ) imode=1;
        else continue;
        // find D* decay products
        const Particle & Dstar = BB.decayProducts()[ix].at(-sign*(423-10*imode))[0];
        if (Dstar.children().size()!=2) continue;
        Particle D0;
        int ipi = imode==0 ?  111 : -sign*211;
        if (Dstar.children()[0].pid()   ==  ipi &&
            Dstar.children()[1].abspid()==-sign*421) {
          D0 = Dstar.children()[1];
        }
        else if (Dstar.children()[1].pid() == ipi &&
                 Dstar.children()[0].abspid()==-sign*421) {
          D0 = Dstar.children()[0];
        }
        else {
          continue;
        }
        // rho decay products
        const Particle& rho = BB.decayProducts()[ix].at( sign*213)[0];
        if (rho.children().size()!=2) continue;
        Particle pip;
        if (rho.children()[1].pid()==111 &&
            rho.children()[0].pid()==sign*211)
          pip = rho.children()[0];
        else if (rho.children()[0].pid()==111 &&
                 rho.children()[1].pid()==sign*211)
          pip = rho.children()[1];
        else {
          continue;
        }
        // boost to B rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix]. mom().betaVec());
        FourMomentum prho   = boost.transform(rho.mom());
        FourMomentum pDstar = boost.transform(Dstar.mom());
        FourMomentum pPi = boost.transform(pip.mom());
        FourMomentum pD  = boost.transform(D0   .mom());
        const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
        pD = boost2.transform(pD);
        double cD = pD.p3().unit().dot(pDstar.p3().unit());
        Vector3 trans1 = pD.p3().unit() - cD*pDstar.p3().unit();
        _h[imode][0]->fill(cD);
        const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
        pPi = boost3.transform(pPi);
        double cPi = pPi.p3().unit().dot(prho.p3().unit());
        Vector3 trans2 = pPi.p3().unit()-cPi*prho.p3().unit();
        _h[imode][1]->fill(cPi);
        double chi = atan2(trans1.cross(trans2).dot(prho.p3().unit()),trans1.dot(trans2));
        _h[imode][2]->fill(chi);
        _p[imode][0]->fill(-0.5*(1-5*sqr(cPi)));
        _p[imode][1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        normalize(_h[ix], 1.0,false);
        Estimate0DPtr tmp;
        book(tmp,1,1,1+ix);
        divide(*_p[ix][0],*_p[ix][1],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    CounterPtr _p[2][2];
    const map<PdgId,unsigned int> mode1   = { {-423,1}, { 213,1}};
    const map<PdgId,unsigned int> mode1CC = { { 423,1}, {-213,1}};
    const map<PdgId,unsigned int> mode2   = { {-413,1}, { 213,1}};
    const map<PdgId,unsigned int> mode2CC = { { 413,1}, {-213,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2003_I611745);

}
