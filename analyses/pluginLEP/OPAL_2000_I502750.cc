// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief rho+/- and omega polarization
  class OPAL_2000_I502750 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2000_I502750);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_ctheta_rho, {0.025, 0.05, 0.1, 0.15, 0.3, 0.6});
      book(_h_ctheta_omega, {0.025, 0.05, 0.1, 0.15, 0.3, 0.6});
      for (size_t i = 0; _h_ctheta_rho->numBins(); ++i) {
        book(_h_ctheta_rho->bin(i+1), "ctheta_rho_"+to_string(i), 20, -1.0, 1.0);
        book(_h_ctheta_omega->bin(i+1), "ctheta_omega_"+to_string(i), 20, -1.0, 1.0);
      }
      book(_h_ctheta_omega_all, "ctheta_omega_all",20,-1.,1.);
    }

    pair<double,double> calcRho(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins() ) {
	double Oi = bin.sumW();
	if(Oi==0.) continue;
	double ai = 0.25*(bin.xMax()*(3.-sqr(bin.xMax())) - bin.xMin()*(3.-sqr(bin.xMin())));
	double bi = 0.75*(bin.xMin()*(1.-sqr(bin.xMin())) - bin.xMax()*(1.-sqr(bin.xMax())));
	double Ei = bin.errW();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    bool findOmegaDecay(Particle omega,Particles & pi0, Particles & pip, Particles & pim) {
      for(const Particle & child : omega.children()) {
	if(child.pid()==211)
	  pip.push_back(child);
	else if(child.pid()==-211)
	  pim.push_back(child);
	else if(child.pid()==111)
	  pi0.push_back(child);
	else if(!child.children().empty()) {
	  if(!findOmegaDecay(child,pi0,pip,pim)) return false;
	}
	else
	  return false;
      }
      return true;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // loop over rho and omega mesons
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==213 || Cuts::abspid==223)) {
	double xE = p.momentum().t()/meanBeamMom;
	Vector3 e1z = p.momentum().p3().unit();
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	if(p.abspid()==213) {
	  if(p.children().size()!=2) continue;
	  int sign = p.pid()/213;
	  Particle pion;
	  if(p.children()[0].pid()==sign*211 && p.children()[1].pid()==111) {
	    pion = p.children()[0];
	  }
	  else if(p.children()[1].pid()==sign*211 && p.children()[0].pid()==111) {
	    pion = p.children()[1];
	  }
	  else
	    continue;
	  Vector3 axis1 = boost.transform(pion.momentum()).p3().unit();
	  double ctheta = e1z.dot(axis1);
	  _h_ctheta_rho->fill(xE,ctheta);
	}
	else {
	  Particles pi0,pip,pim;
	  bool three_pi = findOmegaDecay(p,pi0,pip,pim);
	  if(!three_pi || pi0.size()!=1 || pip.size()!=1 || pim.size()!=1)
	    continue;
	  Vector3 v1 = boost.transform(pi0[0].momentum()).p3().unit();
	  Vector3 v2 = boost.transform(pip[0].momentum()).p3().unit();
	  Vector3 norm = v1.cross(v2).unit();
	  double ctheta = e1z.dot(norm);
	  _h_ctheta_omega->fill(xE,ctheta);
	  if(xE>0.025) _h_ctheta_omega_all->fill(ctheta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      vector<double> x = {0.025,0.05,0.1,0.15,0.3,0.6};
      Estimate1DPtr h_rho  ;
      book(h_rho, 1,1,1);
      Estimate1DPtr h_omega;
      book(h_omega, 2,1,1);
      for (size_t ix=1; ix<_h_ctheta_rho->numBins()+1; ++ix) {
        // rho
        normalize(_h_ctheta_rho->bin(ix));
        pair<double,double> rho00 = calcRho(_h_ctheta_rho->bin(ix));
        h_rho->bin(ix+1).set(rho00.first, rho00.second);
        // omega
        normalize(_h_ctheta_omega->bin(ix));
        rho00 = calcRho(_h_ctheta_omega->bin(ix));
        h_omega->bin(ix+1).set(rho00.first, rho00.second);
      }
      // omega over whole range
      Estimate1DPtr h_omega_all;
      book(h_omega_all,2,2,1);
      normalize(_h_ctheta_omega_all);
      pair<double,double> rho00 = calcRho(_h_ctheta_omega_all);
      h_omega_all->bin(1).set(rho00.first, rho00.second);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_ctheta_rho, _h_ctheta_omega;
    Histo1DPtr _h_ctheta_omega_all;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OPAL_2000_I502750);


}
