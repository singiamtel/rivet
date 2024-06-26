// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> p pbar K*
  class BELLE_2008_I778726 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I778726);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 313);
      BB.addStable(-313);
      BB.addStable( 323);
      BB.addStable(-323);
      declare(BB, "BB");
      // histos
      for(unsigned int iy=0;iy<2;++iy) {
        book(_c[iy],"TMP/c_+"+toString(iy+1));
        for (unsigned int ix=0;ix<3;++ix) {
          book(_h[ix][iy],ix+1,1,1+iy);
        }
      }
    }

    bool isK(int id) const {
      return id==321 || id==311 || id==310 || id==130;
    }

    bool isPi(int id) const {
      return id==211 || id==111;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int sign = 1;
        int imode=0;
        if (BB.decaying()[ix].abspid()==511) {
          _c[0]->fill();
          if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
            sign=1;
          }
          else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
            sign=-1;
          }
          else {
            continue;
          }
        }
        else {
          imode=1;
          _c[1]->fill();
          if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
            sign=1;
          }
          else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
            sign=-1;
          }
          else {
            continue;
          }
        }
        int iK = imode==0 ? 313 : 323;
        const Particle& pp    = BB.decayProducts()[ix].at( sign*2212)[0];
        const Particle& pbar  = BB.decayProducts()[ix].at(-sign*2212)[0];
        const Particle& Kstar = BB.decayProducts()[ix].at( sign*iK  )[0];
        FourMomentum ppbar = pp.mom()+pbar.mom();
        double mppbar = ppbar.mass();
        _h[0][imode]->fill(mppbar);
        if (mppbar<2.85) {
          // boost to B rest frame
          LorentzTransform boost =
            LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix]. mom().betaVec());
          FourMomentum pKstar  = boost.transform(Kstar.mom());
          // proton helicity angle
          FourMomentum pProton = boost.transform(pp.mom());
          ppbar = boost.transform(ppbar);
          LorentzTransform boost1 =
            LorentzTransform::mkFrameTransformFromBeta(ppbar.betaVec());
          pProton = boost1.transform(pProton);
          _h[2][imode]->fill(pKstar.p3().unit().dot(pProton.p3().unit()));
          Particle KK;
          if (isK (Kstar.children()[0].abspid()) &&
             isPi(Kstar.children()[1].abspid())) {
            KK = Kstar.children()[0];
          }
          else if (isK (Kstar.children()[1].abspid()) &&
            isPi(Kstar.children()[0].abspid())) {
            KK = Kstar.children()[1];
          }
          else {
            continue;
          }
          FourMomentum pK = boost.transform(KK.mom());
          const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
          pK = boost3.transform(pK);
          FourMomentum pB = boost3.transform(boost.transform(BB.decaying()[ix].mom()));
          const double cosK = -pB.p3().unit().dot(pK.p3().unit());
          _h[1][imode]->fill(cosK);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<3; ++ix) {
        for (unsigned int iy=0; iy<2; ++iy) {
          if (ix==0) scale    (_h[ix][iy], 1e6/ *_c[iy]);
          else       normalize(_h[ix][iy], 1.0, false);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    CounterPtr _c[2];
    const map<PdgId,unsigned int> mode1   = { { 2212,1},{-2212,1}, { 313,1}};
    const map<PdgId,unsigned int> mode1CC = { { 2212,1},{-2212,1}, {-313,1}};
    const map<PdgId,unsigned int> mode2   = { { 2212,1},{-2212,1}, { 323,1}};
    const map<PdgId,unsigned int> mode2CC = { { 2212,1},{-2212,1}, {-323,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2008_I778726);

}
