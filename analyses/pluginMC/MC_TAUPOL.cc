// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  /// @brief Monte Carlo validation observables for tau polarisation
  class MC_TAUPOL : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_TAUPOL);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 0; // default is ditau
      if ( getOption("MODE") == "DITAU" ) { _mode = 0; }
      if ( getOption("MODE") == "SINGLE" ) { _mode = 1; }
      _lothresh = getOption("LOMASS", _mode? 75*GeV : 86.*GeV);
      _hithresh = getOption("HIMASS", _mode? 85*GeV : 96*GeV);

      // Initialise and register projections
      declare(TauFinder(TauDecay::ANY), "taus");

      // Book histograms
      // tau invmass
      book(_h["m_tau"], "h_m_tau", 50,1.,3.);
      book(_h["m_tautau"], "h_m_tautau", 50,35.,420.);
      book(_h["m_B"], "h_m_B", 50,35.,420.);
      book(_h["B_pT"], "h_B_pT", 100,0.,500.);

      // longnitudal spin correlations
      book(_h["m_tau_vis"], "h_m_tau_vis",50,0.,5.);
      book(_h["m_tautau_vis"], "h_m_tautau_vis",50,0.,100.);
      book(_h["m_pipi_vis"], "h_m_pipi_vis",60,0.,100.);
      book(_h["m_pipi_vis_lm"], "h_m_pipi_vis_lm",60,0.,100.);
      book(_h["m_pipi_vis_hm"], "h_m_pipi_vis_hm",60,0.,100.);

      // transverse spin correlations
      book(_h["acoplanarity_plus"], "h_acoplanarity_plus",30,0.,6.2);
      book(_h["acop_pt40"], "h_acop_pt40",30,0.,6.2);
      book(_h["acop_pt100"], "h_acop_pt100",30,0.,6.2);
      book(_h["acop_pt200"], "h_acop_pt200",30,0.,6.2);
      book(_h["acop_pt1000"], "h_acop_pt1000",30,0.,6.2);
      book(_h["acom_pt40"], "h_acom_pt40",30,0.,6.2);
      book(_h["acom_pt100"], "h_acom_pt100",30,0.,6.2);
      book(_h["acom_pt200"], "h_acom_pt200",30,0.,6.2);
      book(_h["acom_pt1000"], "h_acom_pt1000",30,0.,6.2);
      book(_h["acoplanarity_minus"], "h_acoplanarity_minus",30,0.,6.2);

      // leptonic

      book(_h["2B_xlep"], "h_2B_xlep", 50,0.0,1.0);
      book(_h["2B_xel"], "h_2B_xel", 50,0.0,1.0);
      book(_h["2B_xmu"], "h_2B_xmu", 50,0.0,1.0);
      book(_h["xel_B"], "h_xel_B", 50,0.0,1.0);
      book(_h["xmu_B"], "h_xmu_B", 50,0.0,1.0);

      // hadronic

      book(_h["1B_xpi"], "h_1B_xpi", 30,0.0,1.0);
      book(_h["1B_xpi_lm"], "h_1B_xpi_lm", 30,0.0,1.0);
      book(_h["1B_xpi_hm"], "h_1B_xpi_hm", 30,0.0,1.0);
      book(_h["xpi_B"], "h_xpi_B", 30,0.0,1.0);
      book(_h["2B_xrho"], "h_2B_xrho", 50,0.0,1.0);

      // angles
      book(_h["t_pi"], "h_t_pi", 40, -1.0, 1.0);
      book(_h["t_rho"], "h_t_rho", 40, -1.0, 1.0);
      book(_h["pipi_theta"], "h_pipi_theta", 40, -1.0, 1.0);
      book(_ntaus, "h_n_taus", {0, 1, 2, 3, 4});

    }

    void findDecayProducts(const Particle& mother, size_t& nstable,
			   Particles& ep , Particles& em , Particles& nu_e , Particles& nu_ebar,
			   Particles& mup, Particles& mum, Particles& nu_mu, Particles& nu_mubar,
			   Particles& pip, Particles& pim, Particles& pi0,
			   Particles& Kp , Particles& Km , Particles& K0S,   Particles& K0L,
			   Particles& eta, Particles& gamma) {
      for (const Particle& p : mother.children()) {
        int id = p.pid();
        if (id == PID::KPLUS) {
       	  Kp.push_back(p);
          ++nstable;
        }
        else if (id == PID::KMINUS ) {
          Km.push_back(p);
          ++nstable;
        }
        else if (id == PID::PIPLUS) {
          pip.push_back(p);
          ++nstable;
        }
        else if (id == PID::PIMINUS) {
          pim.push_back(p);
          ++nstable;
        }
        else if (id == PID::EPLUS) {
          ep.push_back(p);
          ++nstable;
        }
        else if (id == PID::EMINUS) {
          em.push_back(p);
          ++nstable;
        }
        else if (id == PID::NU_E) {
          nu_e.push_back(p);
          ++nstable;
        }
        else if (id == PID::NU_EBAR) {
          nu_ebar.push_back(p);
          ++nstable;
        }
        else if (id == PID::NU_MU) {
          nu_mu.push_back(p);
          ++nstable;
        }
        else if (id == PID::NU_MUBAR) {
          nu_mubar.push_back(p);
          ++nstable;
        }
        else if (id == PID::ANTIMUON) {
          mup.push_back(p);
          ++nstable;
        }
        else if (id == PID::MUON) {
          mum.push_back(p);
          ++nstable;
        }
        else if (id == PID::PI0) {
          pi0.push_back(p);
          ++nstable;
        }
        else if (id == PID::K0S) {
          K0S.push_back(p);
          ++nstable;
        }
        else if (id == PID::K0L) {
          K0L.push_back(p);
          ++nstable;
        }
        else if (id == PID::ETA) {
          eta.push_back(p);
          ++nstable;
        }
        else if (id == PID::PHOTON) {
          gamma.push_back(p);
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, ep, em, nu_e, nu_ebar, mup, mum, nu_mu, nu_mubar,
                            pip, pim, pi0, Kp, Km, K0S, K0L, eta, gamma);
        }
        else {
          ++nstable;
        }
      }
    }


    double angleRF(const FourMomentum& p1, const FourMomentum& p2, const FourMomentum& RF) const {
       LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(RF.betaVec());
       return cos(boost(p1).angle(boost(p2)));
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles taus, t_nu, pions;
      for (const auto& tau : apply<TauFinder>(event, "taus").particlesByPt()) {
        //remove any from hadronic decay (only want prompt taus)
        if (tau.parents()[0].isPrompt()) {
          taus += tau;
          _h["m_tau"]->fill(tau.mom().mass()/GeV);
        }
      }
      //veto ditau events without exactly 2 prompt taus
      if (taus.size() != 2 && _mode == 0) vetoEvent;
      //veto single-tau events without exactly 1 prompt tau
      if (taus.size() != 1 && _mode == 1) vetoEvent;

      FourMomentum tautaumom, taunumom, phomom, Bmom;
      if (taus.size() == 2 && _mode == 0) {
	      tautaumom = taus[0].mom() + taus[1].mom();
        phomom = sum(scanAncestors(taus,22), Kin::p4, FourMomentum());
        Bmom = tautaumom + phomom;
        _h["m_tautau"]->fill(tautaumom.mass());
      }
      if (taus.size()==1 && _mode == 1) {
        t_nu = scanAncestors(taus,16);
        phomom = sum(scanAncestors(taus,22), Kin::p4, FourMomentum());
        if(t_nu.size() != 1)vetoEvent;
        if(taus[0].pid()*t_nu[0].pid() > 0) vetoEvent;
        taunumom = taus[0].mom() + t_nu[0].mom();
        Bmom = taunumom + phomom;
      }
      _ntaus->fill(int(taus.size()));
      _h["m_B"]->fill(Bmom.mass());
      _h["B_pT"]->fill(Bmom.pT());
       // Check opposite charges

      // Calculate invariant mass of tautau
      FourMomentum tt_mom_vis;
      for (const auto& t : taus) {
        FourMomentum t_mom_vis;
        for (const auto& c : t.children()) {
          if (c.isVisible()) {
            tt_mom_vis += c.mom();
            t_mom_vis += c.mom();
          }
        }
        _h["m_tau_vis"]->fill(t_mom_vis.mass());
      }

      _h["m_tautau_vis"]->fill(tt_mom_vis.mass());
      size_t OS = 2;
      if (Bmom.mass() < _lothresh)  OS = 0;
      else if (Bmom.mass() < _hithresh)  OS = 1;

      for (const Particle& tau :taus) {
        size_t nstable(0);
        Particles ep,em,nu_e,nu_ebar,mup,mum,nu_mu,nu_mubar;
        Particles pip, pim, pi0, Kp , Km, K0S, K0L, eta,gamma;
        findDecayProducts(tau, nstable,ep,em,nu_e,nu_ebar,mup,mum,nu_mu,nu_mubar,
                          pip, pim, pi0, Kp, Km, K0S, K0L,eta,gamma);
        if (tau.pid() < 0) {
          swap(pim,pip);
          swap(Kp,Km);
          swap(em,ep);
          swap(mum,mup);
          swap(nu_e ,nu_ebar );
          swap(nu_mu,nu_mubar);
        }
        // classify by number of stable decay products
        FourMomentum tmom = tau.mom();
        const LorentzTransform tboost = LorentzTransform::mkFrameTransformFromBeta(tmom.betaVec());
        if (nstable==2) {
          if (pim.size()==1) {
            const double xpi = pim[0].mom().p()/tau.mom().p();
            const FourMomentum pimom = pim[0].mom();
            const double angle = cos(tboost(pimom).angle(tmom));
            _h["xpi_B"]->fill(pim[0].mom().p()/Bmom.p());
            _h["t_pi"]->fill(angle);
            if (OS==0)      _h["1B_xpi_lm"]->fill(xpi);
            else if (OS==1) _h["1B_xpi"]->fill(xpi);
            else            _h["1B_xpi_hm"]->fill(xpi);
            pions.push_back(pim[0]);
          }
        }
        else if (nstable==3 ) {
          if (em.size()==1 && nu_ebar.size()==1) {
            const double xlep = em[0].mom().p()/tau.mom().p();
            _h["2B_xlep"]->fill(xlep);
            _h["2B_xel"]->fill(xlep);
            _h["xel_B"]->fill(em[0].mom().p()/Bmom.p());
          }
          else if (mum.size()==1 && nu_mubar.size()==1) {
            const double xlep = mum[0].mom().p()/tau.mom().p();
            _h["2B_xlep"]->fill(xlep);
            _h["2B_xmu"]->fill(xlep);
            _h["xmu_B"]->fill(mum[0].mom().p()/Bmom.p());
          }
          else if (pim.size()==1 && pi0.size()==1) {
            const FourMomentum ptot = pim[0].mom()+pi0[0].mom();
            const double xrho = ptot.p()/tau.mom().p();
            _h["2B_xrho"]->fill(xrho);
            const double angle = cos(tboost(ptot).angle(tmom));
            _h["t_rho"]->fill(angle);
          }
        } //end nstable = 3

      } //end tau loop

      // TAU LONGNITUDAL SPIN EFFECTS PART
      //fill tautau->pipi inv mass
      if (pions.size() == 2){
        FourMomentum pi1mom = pions[0].mom();
        FourMomentum pi2mom = pions[1].mom();
        FourMomentum pipi = pions[0].mom() + pions[1].mom();
        if (OS==0)      _h["m_pipi_vis_lm"]->fill(pipi.mass()/GeV);
        else if (OS==1) _h["m_pipi_vis"]->fill(pipi.mass()/GeV);
        else            _h["m_pipi_vis_hm"]->fill(pipi.mass()/GeV);
        _h["pipi_theta"]->fill(angleRF(pi1mom,pi2mom,Bmom));
      }

  		// TAU TRANSVERSE SPIN EFFECTS PART
      if (pions.size() == 2 && _mode != 1 && OS == 1){

        FourMomentum tautau = taus[0].mom() + taus[1].mom();

  	    // calculate angles between (tau, tau_pi) and (tau2, tau2_pi) planes
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(tautau.betaVec());
	      taus[0] = Particle(PID::ANY, boost.transform(taus[0].mom()));
	      pions[0] = Particle(PID::ANY, boost.transform(pions[0].mom()));
        taus[1] = Particle(PID::ANY, boost.transform(taus[1].mom()));
	      pions[1] = Particle(PID::ANY, boost.transform(pions[1].mom()));

        // calculate the angle
        Vector3 t0xp0 = taus[0].p3().cross(pions[0].p3()).unit();
        Vector3 t1xp1 = taus[1].p3().cross(pions[1].p3()).unit();
        double phi_star_sign = acos(t0xp0.dot(t1xp1));

        // extend acoplanarity angle to range (0, 2pi)
        Vector3 p0xp1 = pions[0].p3().cross(pions[1].p3()).unit();
        double sign = p0xp1.dot(taus[0].p3().unit());
        if (sign > 0)  phi_star_sign = 2* M_PI - phi_star_sign;


        // alpha angles in the lab frame
        Vector3  tlv3_Ez(0.0, 0.0, 1.0 );
        Vector3 cross1_pi = tlv3_Ez.cross(taus[0].p3()).unit();
        Vector3 cross2_pi = pions[0].p3().cross(taus[0].p3()).unit();
        const double alpha_pi = acos(abs(cross1_pi.dot(cross2_pi)));

        if (alpha_pi > M_PI/4.) {
          _h["acoplanarity_plus"]->fill(phi_star_sign);
          if (Bmom.pT() <  40*GeV)       _h["acop_pt40"]->fill(phi_star_sign);
          else if (Bmom.pT() < 100*GeV)  _h["acop_pt100"]->fill(phi_star_sign);
          else if (Bmom.pT() < 200*GeV)  _h["acop_pt200"]->fill(phi_star_sign);
          else                           _h["acop_pt1000"]->fill(phi_star_sign);
        }

        if (alpha_pi < M_PI/4.) {
          _h["acoplanarity_minus"]->fill(phi_star_sign);
          if (Bmom.pT() <  40*GeV)       _h["acom_pt40"]->fill(phi_star_sign);
          else if (Bmom.pT() < 100*GeV)  _h["acom_pt100"]->fill(phi_star_sign);
          else if (Bmom.pT() < 200*GeV)  _h["acom_pt200"]->fill(phi_star_sign);
          else                           _h["acom_pt1000"]->fill(phi_star_sign);
        }

      } // end of tau spin corr

    } //end void function


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
      normalize(_ntaus);
    }

    //@}

    Particles scanAncestors(const Particles& taus, int id){
      Particles rtn;
      for (const Particle& tau : taus) {
        Particle thistau = tau;
        while (true) {
          // check the tau has parents first else end loop
          if (thistau.parents().empty())  break;
          //loop over taus siblings
          for (const Particle& sibling : thistau.parents()[0].children()) {
            // find ID
            if (sibling.abspid() == id) {
              // add particles to rtn vector
              bool duplicate = true;
              for (const auto& p : rtn) {
                if (sibling.isSame(p)) duplicate = false;
              }
              if (duplicate) rtn += sibling;
            }
          }
          // want to loop over taus parent's siblings next
          thistau = thistau.parents()[0];
          // if taus parent is not a tau end the loop
          if (thistau.abspid() != PID::TAU)  break;
        }
      }
      return rtn;
    }

    /// @name Histograms
    //@{

    // histograms for leptonic decay
    map<string,Histo1DPtr> _h;

    BinnedHistoPtr<int> _ntaus;

    size_t _mode;

    double _lothresh, _hithresh;

    //@}

  };


  RIVET_DECLARE_PLUGIN(MC_TAUPOL);
}
