// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B-> p Lambdabar X decays
  class BELLE_2007_I748975 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I748975);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PI0);
      BB.addStable( 3122);
      BB.addStable(-3122);
      declare(BB, "BB");
      // histos
      for (unsigned int iy=0; iy<3; ++iy) {
        book(_h_mass [iy],1,1,1+iy);
        book(_h_cosp [iy],2,1,1+iy);
        book(_h_alpha[iy][0],3,1,1+iy);
        book(_h_alpha[iy][1],"TMP/h_alpha_"+toString(1+iy));
        book(_h_E    [iy][0],4,1,1+iy);
        book(_h_E    [iy][1],"TMP/h_E_"+toString(1+iy));
      }
      for (unsigned int iy=0;iy<2;++iy) {
        book(_c[iy],"TMP/c_"+toString(iy+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int imode=0, sign=1, idOther=22;
        if (BB.decaying()[ix].abspid()==521) {
          _c[0]->fill();
          if (BB.modeMatches(ix,3,mode1)) {
            imode= 0;
            sign = 1;
          }
          else if (BB.modeMatches(ix,3,mode1CC)) {
            imode= 0;
            sign =-1;
          }
          else if (BB.modeMatches(ix,3,mode2)) {
            imode= 1;
            sign = 1;
            idOther= 111;
          }
          else if (BB.modeMatches(ix,3,mode2CC)) {
            imode= 1;
            sign =-1;
            idOther= 111;
          }
          else continue;
        }
        else {
          _c[1]->fill();
          if (BB.modeMatches(ix,3,mode3)) {
            imode= 2;
            sign = 1;
            idOther=-211;
          }
          else if (BB.modeMatches(ix,3,mode3CC)) {
            imode= 2;
            sign =-1;
            idOther= 211;
          }
          else continue;
        }
        const Particle& pp = BB.decayProducts()[ix].at( sign*2212)[0];
        const Particle& lbar  = BB.decayProducts()[ix].at(-sign*3122)[0];
        const Particle& other = BB.decayProducts()[ix].at(idOther)[0];
        FourMomentum pbaryon = pp.mom()+lbar.mom();
        double mass = pbaryon.mass();
        _h_mass[imode]->fill(mass);
        // boost to B rest frame
        LorentzTransform boost =
          LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix]. mom().betaVec());
        FourMomentum pLam = boost.transform(lbar.mom());
        pbaryon = boost.transform(pbaryon);
        FourMomentum pProton = boost.transform(pp.mom());
        _h_E[imode][0]->fill(pLam.E());
        _h_E[imode][1]->fill();
        if (mass<2.8) {
          Vector3 axis = boost.transform(other.mom()).p3().unit();
          LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(pbaryon.betaVec());
          double cosp = axis.dot(boost1.transform(pProton).p3().unit());
          _h_cosp[imode]->fill(cosp);
        }
        // finally the alpha parameter
        // Lambda decay products
        if (lbar.children().size()!=2) continue;
        Particle pbar;
        if (lbar.children()[0].pid()==-2212 &&
           lbar.children()[1].pid()== 211) {
          pbar = lbar.children()[0];
        }
        else if (lbar.children()[1].pid()==-sign*2212 &&
          lbar.children()[0].pid()== sign*211) {
          pbar = lbar.children()[1];
        }
        else {
          continue;
        }
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pLam.betaVec());
        Vector3 axis1 = - boost2.transform(boost.transform(BB.decaying()[ix].mom())).p3().unit();
        Vector3 axis2 = - boost2.transform(boost.transform(pbar.mom())).p3().unit();
        double cTheta = axis1.dot(axis2);
        _h_alpha[imode][0]->fill(1.5*cTheta);
        _h_alpha[imode][1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_cosp, 1.0, false);
      scale(_h_mass[0],1e6/ *_c[0]);
      scale(_h_mass[1],1e6/ *_c[0]);
      scale(_h_mass[2],1e6/ *_c[1]);
      for(unsigned int ix=0;ix<3;++ix) {
        scale(_h_alpha[ix][0], 1./ *_h_alpha[ix][1]);
        scale(_h_E    [ix][0], 1./ *_h_E    [ix][1]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3],_h_cosp[3];
    CounterPtr _h_alpha[3][2],_h_E[3][2];
    CounterPtr _c[2];
    const map<PdgId,unsigned int> mode1   = { { 2212,1}, {-3122,1}, { 22,1} };
    const map<PdgId,unsigned int> mode1CC = { {-2212,1}, { 3122,1}, { 22,1} };
    const map<PdgId,unsigned int> mode2   = { { 2212,1}, {-3122,1}, { 111,1} };
    const map<PdgId,unsigned int> mode2CC = { {-2212,1}, { 3122,1}, { 111,1} };
    const map<PdgId,unsigned int> mode3   = { { 2212,1}, {-3122,1}, {-211,1} };
    const map<PdgId,unsigned int> mode3CC = { {-2212,1}, { 3122,1}, { 211,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I748975);

}
