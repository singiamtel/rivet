// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> omega pi/K
  class BELLE_2006_I735859 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I735859);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      BB.addStable(223);
      declare(BB, "BB");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }

    void findChildren(const Particle& p, Particles& pim, Particles& pip,
                      Particles& pi0, unsigned int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.pid()==PID::PIPLUS) {
          pip.push_back(child);
          ncount+=1;
        }
        else if (child.pid()==PID::PIMINUS) {
          pim.push_back(child);
          ncount+=1;
        }
        else if (child.pid()==PID::PI0) {
          pi0.push_back(child);
          ncount+=1;
        }
        else if (child.children().empty()) {
          ncount+=1;
        }
        else {
          findChildren(child,pim,pip,pi0,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        unsigned int imode=0;
        if (BB.modeMatches(ix,2,mode1) || BB.modeMatches(ix,2,mode1CC)) {
          imode=0;
        }
        else if (BB.modeMatches(ix,2,mode2) || BB.modeMatches(ix,2,mode2CC)) {
          imode=1;
        }
        else if (BB.modeMatches(ix,2,mode3) || BB.modeMatches(ix,2,mode3CC)) {
          imode=2;
        }
        else {
          continue;
        }
        const Particle& omega = BB.decayProducts()[ix].at(223)[0];
        // children of the omega
        unsigned int ncount=0;
        Particles pip,pim,pi0;
        findChildren(omega,pim,pip,pi0,ncount);
        if (ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) continue;
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].mom().betaVec());
        FourMomentum pomega = boost1.transform(omega.mom());
        FourMomentum pPip   = boost1.transform(pip[0].mom());
        FourMomentum pPim   = boost1.transform(pim[0].mom());
        // boost to omega frame
        LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pomega.betaVec());
        pPip   = boost3.transform(pPip);
        pPim   = boost3.transform(pPim);
        Vector3 axisOmega = pPip.p3().cross(pPim.p3()).unit();
        // helicity angle
        _h[imode]->fill(abs(pomega.p3().unit().dot(axisOmega)));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    const map<PdgId,unsigned int> mode1   = { { 223,1}, { 321,1}};
    const map<PdgId,unsigned int> mode1CC = { { 223,1}, {-321,1}};
    const map<PdgId,unsigned int> mode2   = { { 223,1}, { 211,1}};
    const map<PdgId,unsigned int> mode2CC = { { 223,1}, {-211,1}};
    const map<PdgId,unsigned int> mode3   = { { 223,1}, { 130,1}};
    const map<PdgId,unsigned int> mode3CC = { { 223,1}, { 310,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2006_I735859);

}
