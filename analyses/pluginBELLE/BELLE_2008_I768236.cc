// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D**  l nu
  class BELLE_2008_I768236 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I768236);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(411); BB.addStable(-411);
      BB.addStable(421); BB.addStable(-421);
      BB.addStable(413); BB.addStable(-413);
      BB.addStable(423); BB.addStable(-423);
      BB.addStable(PID::PI0);
      declare(BB, "BB");
      // hists
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h_mass[ix], 1, 1, 1+ix);
        if (ix==3 || ix==1) continue;
        book(_h_hel [ix], 2, 1, 1+ix);
        book(_h_w   [ix], 3, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> mode1[4] = {{ {-411,1}, { 211,1}, {-11,1}, { 12,1}},
                                                       { { 411,1}, {-211,1}, { 11,1}, {-12,1}},
                                                       { {-411,1}, { 211,1}, {-13,1}, { 14,1}},
                                                       { { 411,1}, {-211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode2[4] = {{ {-413,1}, { 211,1}, {-11,1}, { 12,1}},
                                                       { { 413,1}, {-211,1}, { 11,1}, {-12,1}},
                                                       { {-413,1}, { 211,1}, {-13,1}, { 14,1}},
                                                       { { 413,1}, {-211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode3[4] = {{ {-421,1}, {-211,1}, {-11,1}, { 12,1}},
                                                       { { 421,1}, { 211,1}, { 11,1}, {-12,1}},
                                                       { {-421,1}, {-211,1}, {-13,1}, { 14,1}},
                                                       { { 421,1}, { 211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode4[4] = {{ {-423,1}, {-211,1}, {-11,1}, { 12,1}},
                                                       { { 423,1}, { 211,1}, { 11,1}, {-12,1}},
                                                       { {-423,1}, {-211,1}, {-13,1}, { 14,1}},
                                                       { { 423,1}, { 211,1}, { 13,1}, {-14,1}}};
      // loop over B mesons
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
      	for (unsigned int il=0; il<4; ++il) {
          int iD,imode=0;
          if ( BB.modeMatches(ix,4,mode1[il]) ) {
            imode=0;
            iD=411;
          }
          else if ( BB.modeMatches(ix,4,mode2[il]) ) {
            imode=1;
            iD=413;
          }
          else if ( BB.modeMatches(ix,4,mode3[il]) ) {
            imode=2;
            iD=421;
          }
          else if ( BB.modeMatches(ix,4,mode4[il]) ) {
            imode=3;
            iD=423;
          }
          else continue;
          int sign = il%2==0 ? 1 : -1;
          int ipi = -sign*211;
          iD *= -sign;
          if(imode<2) ipi*=-1;
          const Particle & pi1= BB.decayProducts()[ix].at( ipi)[0];
          const Particle & DD = BB.decayProducts()[ix].at( iD )[0];
          FourMomentum pHad = pi1.mom()+DD.mom();
          double mass = pHad.mass();
          _h_mass[imode]->fill(mass);
          // just the D pi modes for helicty and w
          if (imode%2!=0) continue;
          // check ingtermediate
          int inter=-1;
          for (unsigned int iz=0;iz<BB.decaying()[ix].children().size();++iz) {
            if (BB.decaying()[ix].children()[iz].abspid()==10411 ||
                BB.decaying()[ix].children()[iz].abspid()==10421) {
              inter=0;
              break;
            }
            else if (BB.decaying()[ix].children()[iz].abspid()==415 ||
                     BB.decaying()[ix].children()[iz].abspid()==425) {
              inter=2;
              break;
            }
          }
          if (inter<0) continue;
          // First boost all relevant momenta into the B-rest frame
          const LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].mom().betaVec());
          // Momenta in B rest frame:
          FourMomentum pDSS = boost1.transform(pHad);
          FourMomentum pW   = boost1.transform(BB.decaying()[ix].mom() - pHad);
          FourMomentum ppi  = boost1.transform(pi1);
          const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDSS.betaVec());
          double ctheta = pW.p3().unit().dot(boost2.transform(ppi).p3().unit());
          _h_hel[inter]->fill(abs(ctheta));
          double mB2 = BB.decaying()[ix].mom().mass2();
          double mD2 = sqr(mass);
          double w = (mB2 + mD2 - pW.mass2())/ (2. * sqrt(mB2) * sqrt(mD2) );
          _h_w[inter]->fill(w);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass, 1.0, false);
      normalize(_h_hel, 1.0, false);
      normalize(_h_w, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[4],_h_hel[3],_h_w[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2008_I768236);

}
