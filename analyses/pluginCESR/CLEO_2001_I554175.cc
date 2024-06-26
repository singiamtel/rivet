// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D(*) pi+pi-pi-pi0
  class CLEO_2001_I554175 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2001_I554175);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PI0);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( 423);
      BB.addStable(-423);
      BB.addStable( 411);
      BB.addStable(-411);
      BB.addStable( 421);
      BB.addStable(-421);
      declare(BB, "BB");
      // // histos
      for (unsigned int ix=0; ix<6; ++ix) {
      	book(_h[ix],1+ix,1,1);
      }
      for (unsigned int ix=0; ix<3; ++ix) {
      	book(_h_angle[ix],7,1,1+ix);
      }
      book(_h_sum,8,1,1);
      for (unsigned int ix=0;ix<2;++ix) {
      	book(_c[ix],"TMP/nB_"+toString(ix+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over particles
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      int imode = -1;
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
        if (BB.decaying()[ix].abspid()==511) _c[0]->fill();
        else                                 _c[1]->fill();
        if ((sign== 1 && BB.modeMatches(ix,5,mode1) ) ||
            (sign==-1 && BB.modeMatches(ix,5,mode1CC) )) {
          imode=0;
        }
        else if ((sign== 1 && BB.modeMatches(ix,5,mode2) ) ||
                 (sign==-1 && BB.modeMatches(ix,5,mode2CC) )) {
          imode=1;
        }
        else if ((sign== 1 && BB.modeMatches(ix,5,mode3) ) ||
                 (sign==-1 && BB.modeMatches(ix,5,mode3CC)))
          imode=2;
        else if ((sign== 1 && BB.modeMatches(ix,5,mode4) ) ||
                 (sign==-1 && BB.modeMatches(ix,5,mode4CC) )) {
          imode=3;
        }
        else {
          continue;
        }
        const Particles& pip = BB.decayProducts()[ix].at( sign*211);
        const Particle & pim = BB.decayProducts()[ix].at(-sign*211)[0];
        const Particle & pi0 = BB.decayProducts()[ix].at(      111)[0];
        FourMomentum pOmegaPi = pip[0].mom()+pip[1].mom()+pim.mom()+pi0.mom();
        const double mHad = pOmegaPi.mass();
        if (imode==0) {
          _h[0]->fill(mHad);
        }
        // find the children of the omega
        Particles omegaDec;
        for (const Particle& p : {pip[0],pip[1],pim,pi0} ) {
          Particle parent = p;
          while (parent.parents()[0].pid()!=BB.decaying()[ix].pid()) {
            parent = parent.parents()[0];
            if (parent.pid()==223) {
              omegaDec.push_back(p);
              break;
            }
          }
        }
        if (omegaDec.size()!=3) continue;
        if (imode<2)  _h[imode+1]->fill(mHad);
        else          _h[5      ]->fill(mHad);
        if (imode==1) continue;
        _h_sum->fill(mHad);
        // boost to B rest frame
        LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].mom().betaVec());
        /// D star angles
        if (imode==0) {
          const Particle& Dstar = BB.decayProducts()[ix].at(-sign*413)[0];
          if (Dstar.children().size()==2) {
            Particle D0;
            if (Dstar.children()[0].pid()    ==-sign*211 &&
                Dstar.children()[1].abspid() ==-sign*421) {
              D0 = Dstar.children()[1];
            }
            else if (Dstar.children()[1].pid() ==-sign*211 &&
                     Dstar.children()[0].abspid() ==-sign*421) {
              D0 = Dstar.children()[0];
            }
            // if right decay mode
            if (D0.abspid()==421) {
              FourMomentum pDstar = boostB.transform(Dstar.mom());
              FourMomentum pD0    = boostB.transform(D0   .mom());
              LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
              pD0 = boost2.transform(pD0);
              double c1 = pD0.p3().unit().dot(pDstar.p3().unit());
              _h[3]->fill(c1);
            }
          }
        }
        // omega momenta in B rest frame
        FourMomentum pOmega;
        for (const Particle & p : omegaDec) pOmega+=p.mom();
        pOmega = boostB.transform(pOmega);
        // boost to A rest frame
        pOmegaPi = boostB.transform(pOmegaPi);
        FourMomentum pDstar = boostB.transform(BB.decaying()[ix].mom()) - pOmegaPi;
        LorentzTransform boostWpi = LorentzTransform::mkFrameTransformFromBeta(pOmegaPi.betaVec());
        pOmega = boostWpi.transform(pOmega);
        Vector3 axisW   = pOmega.p3().unit();
        Vector3 axisWpi = pOmegaPi.p3().unit();
        double cA = axisW.dot(axisWpi);
        // omega angles
        LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
        pOmegaPi = boostW.transform(boostWpi.transform(pOmegaPi));
        pDstar = boostW.transform(boostWpi.transform(pDstar));
        axisWpi = pOmegaPi.p3().unit();
        FourMomentum ppip = boostW.transform(boostWpi.transform(boostB.transform(omegaDec[0])));
        FourMomentum ppim = boostW.transform(boostWpi.transform(boostB.transform(omegaDec[1])));
        Vector3 nW = ppip.p3().cross(ppim.p3()).unit();
        double cTheta1 = axisWpi.dot(nW);
        if (imode==0) {
          _h[4]->fill(cTheta1);
          continue;
        }
        _h_angle[0]->fill(cA);
        _h_angle[1]->fill(cTheta1);
        Vector3 transW = nW-cTheta1*axisWpi;
        Vector3 transD = pDstar.p3().unit()-pDstar.p3().unit().dot(axisWpi)*axisWpi;
        double phi = abs(atan(transW.cross(transD).dot(axisWpi)/ transW.dot(transD)));
        _h_angle[2]->fill(phi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // first hist is differential BR
      scale(_h[0], 1./ *_c[0]);
      // rest are unit normalized
      for(unsigned int ix=1;ix<6;++ix)
        normalize(_h[ix], 1.0, false);
      normalize(_h_angle, 1.0, false);
      normalize(_h_sum, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6],_h_angle[3],_h_sum;
    CounterPtr _c[2];
    const map<PdgId,unsigned int> mode1   = { { -413,1}, { 211,2}, {-211,1}, {111,1}};
    const map<PdgId,unsigned int> mode1CC = { {  413,1}, {-211,2}, { 211,1}, {111,1}};
    const map<PdgId,unsigned int> mode2   = { { -423,1}, { 211,2}, {-211,1}, {111,1}};
    const map<PdgId,unsigned int> mode2CC = { {  423,1}, {-211,2}, { 211,1}, {111,1}};
    const map<PdgId,unsigned int> mode3   = { { -411,1}, { 211,2}, {-211,1}, {111,1}};
    const map<PdgId,unsigned int> mode3CC = { {  411,1}, {-211,2}, { 211,1}, {111,1}};
    const map<PdgId,unsigned int> mode4   = { { -421,1}, { 211,2}, {-211,1}, {111,1}};
    const map<PdgId,unsigned int> mode4CC = { {  421,1}, {-211,2}, { 211,1}, {111,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2001_I554175);

}
