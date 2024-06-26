// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi+- K+- K0S pbar spectra continuum and Upsilon 1s
  class ARGUS_1989_I276860 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1989_I276860);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_pi_ups1_p,  5, 1, 1);
      book(_h_pi_cont_p,  5, 1, 2);
      book(_h_pi_ups1_z,  9, 1, 1);
      book(_h_pi_cont_z,  9, 1, 2);

      book(_h_Kp_ups1_p,  6, 1, 1);
      book(_h_Kp_cont_p,  6, 1, 2);
      book(_h_Kp_ups1_z, 10, 1, 1);
      book(_h_Kp_cont_z, 10, 1, 2);

      book(_h_KS_ups1_p,  7, 1, 1);
      book(_h_KS_cont_p,  7, 1, 2);
      book(_h_KS_ups1_z, 11, 1, 1);
      book(_h_KS_cont_z, 11, 1, 2);

      book(_h_pt_ups1_p,  8, 1, 1);
      book(_h_pt_cont_p,  8, 1, 2);
      book(_h_pt_ups1_z, 12, 1, 1);
      book(_h_pt_cont_z, 12, 1, 2);

      // Counters
      book(_n_PiA,1,1,1);
      book(_n_PiB,1,1,2);
      book(_n_KS ,2,1,1);
      book(_n_Kp ,3,1,1);
      book(_n_ptA,4,1,1);
      book(_n_ptB,4,1,2);
      
      // weights
      book(_weightSum_cont,"/TMP/weightSum_cont");
      book(_weightSum_Ups1,"/TMP/weightSum_Ups1");
    }

   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
	if(id == 211 || id == 2212 || id == 310|| id==130 || id == 321) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      const FinalState & fs = apply<FinalState>(event, "FS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      // Continuum
      if (upsilons.empty()) {
        MSG_DEBUG("No Upsilons found => continuum event");
        _weightSum_cont->fill();
        // stable particles
        for (const Particle& p : fs.particles(Cuts::abspid==211 or Cuts::abspid==2212 or Cuts::abspid==321 )) {
          const int id = abs(p.pid());
          const double xE = 2.*p.E()/sqrtS();
          const double modp = p.p3().mod();
          const double beta = modp / p.E();
          int idMom = !p.parents().empty() ? abs(p.parents()[0].pid()) : 0;
          if(id==211) {
            // not from K0S or Lambda decays
            if(idMom!=310 && idMom != 3122) {
              _h_pi_cont_p->fill(modp);
              _h_pi_cont_z->fill(xE  ,1./beta);
              _n_PiA->fill("9.98"s);
            }
            _n_PiB->fill("9.98"s);
          }
          else if(id==321) {
            _h_Kp_cont_p->fill(modp);
            _h_Kp_cont_z->fill(xE  ,1./beta);
            _n_Kp->fill("9.98"s);
          }
          else if(id==2212) {
            // not from K0S or Lambda decays
            if(idMom!=310 && idMom != 3122) {
              _h_pt_cont_p->fill(modp);
              _h_pt_cont_z->fill(xE  ,1./beta);
              _n_ptA->fill("9.98"s);
            }
            _n_ptB->fill("9.98"s);
          }
        }
        // Unstable particles
        for (const Particle& p : ufs.particles(Cuts::pid==310 || Cuts::pid==130)) {
          const double xE = 2.*p.E()/sqrtS();
          const double modp = p.p3().mod();
          const double beta = modp / p.E();
          _h_KS_cont_p->fill(modp);
          _h_KS_cont_z->fill(xE  ,1./beta);
          _n_KS->fill("9.98"s);
        }
      }
      // Upsilon(s) found
      else {
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
          _weightSum_Ups1->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups, unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
          // loop over decay products
          for(const Particle& p : unstable) {
            const int id = abs(p.pid());
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xE = 2.*p2.E()/mass;
            const double modp = p2.p3().mod();
            const double beta = modp / p2.E();
            int idMom = !p.parents().empty() ? abs(p.parents()[0].pid()) : 0;
            if(id==211) {
              // not from K0S or Lambda decays
              if(idMom!=310 && idMom != 3122) {
        	_h_pi_ups1_p->fill(modp);
        	_h_pi_ups1_z->fill(xE  ,1./beta);
        	_n_PiA->fill("9.46"s);
              }
              _n_PiB->fill("9.46"s);
            }
            else if(id==321) {
              _h_Kp_ups1_p->fill(modp);
              _h_Kp_ups1_z->fill(xE  ,1./beta);
              _n_Kp->fill("9.46"s);
            }
            else if(id==2212) {
              // not from K0S or Lambda decays
              if(idMom!=310 && idMom != 3122) {
        	_h_pt_ups1_p->fill(modp);
        	_h_pt_ups1_z->fill(xE  ,1./beta);
        	_n_ptA->fill("9.46"s);
              }
              _n_ptB->fill("9.46"s);
            }
            else if(id==310 || id==130) {
              _h_KS_ups1_p->fill(modp);
              _h_KS_ups1_z->fill(xE  ,1./beta);
              _n_KS->fill("9.46"s);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Scale histos
      if (_weightSum_cont->val() > 0.) {
	scale(_h_pi_cont_p, 1./ *_weightSum_cont);
	scale(_h_Kp_cont_p, 1./ *_weightSum_cont);
	scale(_h_KS_cont_p, 1./ *_weightSum_cont);
	scale(_h_pt_cont_p, 1./ *_weightSum_cont);
	scale(_h_pi_cont_z, 1./ *_weightSum_cont);
	scale(_h_Kp_cont_z, 1./ *_weightSum_cont);
	scale(_h_KS_cont_z, 1./ *_weightSum_cont);
	scale(_h_pt_cont_z, 1./ *_weightSum_cont);

      }
      if (_weightSum_Ups1->val() > 0.) {
	scale(_h_pi_ups1_p, 1./ *_weightSum_Ups1);
	scale(_h_Kp_ups1_p, 1./ *_weightSum_Ups1);
	scale(_h_KS_ups1_p, 1./ *_weightSum_Ups1);
	scale(_h_pt_ups1_p, 1./ *_weightSum_Ups1);
	scale(_h_pi_ups1_z, 1./ *_weightSum_Ups1);
	scale(_h_Kp_ups1_z, 1./ *_weightSum_Ups1);
	scale(_h_KS_ups1_z, 1./ *_weightSum_Ups1);
	scale(_h_pt_ups1_z, 1./ *_weightSum_Ups1);
      }
      // Counters
      vector<CounterPtr> scales = {_weightSum_Ups1,_weightSum_cont};
      vector<string> labels={"9.46","9.98"};
      for (unsigned int ix=0;ix<2;++ix) {
        _n_PiA->binAt(labels[ix]).scaleW(1./ scales[ix]->val());
        _n_PiB->binAt(labels[ix]).scaleW(1./ scales[ix]->val());
        _n_Kp->binAt(labels[ix] ).scaleW(1./ scales[ix]->val());
        _n_KS->binAt(labels[ix] ).scaleW(1./ scales[ix]->val());
        _n_ptA->binAt(labels[ix]).scaleW(1./ scales[ix]->val());
        _n_ptB->binAt(labels[ix]).scaleW(1./ scales[ix]->val());
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi_ups1_p, _h_pi_cont_p, _h_pi_ups1_z, _h_pi_cont_z;
    Histo1DPtr _h_Kp_ups1_p, _h_Kp_cont_p, _h_Kp_ups1_z, _h_Kp_cont_z;
    Histo1DPtr _h_KS_ups1_p, _h_KS_cont_p, _h_KS_ups1_z, _h_KS_cont_z;
    Histo1DPtr _h_pt_ups1_p, _h_pt_cont_p, _h_pt_ups1_z, _h_pt_cont_z;

      // Counters
    BinnedHistoPtr<string> _n_PiA, _n_PiB, _n_Kp, _n_KS, _n_ptA,_n_ptB;
    CounterPtr _weightSum_cont,_weightSum_Ups1;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ARGUS_1989_I276860);

}
