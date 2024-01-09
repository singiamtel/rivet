// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /// @brief CDF Run II underlying event in Drell-Yan
  ///
  /// @author Hendrik Hoeth
  ///
  /// Measurement of the underlying event in Drell-Yan
  /// \f$ Z/\gamma^* \to e^+ e^- \f$ and
  /// \f$ Z/\gamma^* \to \mu^+ \mu^- \f$ events. The reconstructed
  /// Z defines the \f$ \phi \f$ orientation. A Z mass window cut is applied.
  ///
  /// @par Run conditions
  ///
  /// @arg \f$ \sqrt{s} = \f$ 1960 GeV
  /// @arg produce Drell-Yan events
  /// @arg Set particles with c*tau > 10 mm stable
  /// @arg Z decay mode: Z -> e+e- and Z -> mu+mu-
  /// @arg gamma decay mode: gamma -> e+e- and gamma -> mu+mu-
  /// @arg minimum invariant mass of the fermion pair coming from the Z/gamma: 70 GeV
  class CDF_2010_I849042 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2010_I849042);


    /// @name Analysis methods
    /// @{

    void init() {

      _mode = 0;
      if ( getOption("MODE") == "DY" ) _mode = 1;
      else if ( getOption("MODE") == "QCD" ) _mode = 2;

      // Set up projections
      const ChargedFinalState cfs(Cuts::abseta < 1.0 && Cuts::pT >= 0.5*GeV);
      const ChargedFinalState clfs(Cuts::abseta < 1.0 && Cuts::pT >= 20*GeV);
      declare(cfs, "CFS");
      declare(ChargedLeptons(clfs), "CL");

      // Final state for the jet finding
      const FinalState fsj(Cuts::abseta < 4.0);
      declare(fsj, "FSJ");
      declare(FastJets(fsj, JetAlg::CDFMIDPOINT, 0.7), "MidpointJets");

      // Book histograms
      if (_mode == 0 || _mode == 1) {
        book(_p["tnchg"],                1, 1, 1);
        book(_p["pnchg"],                1, 1, 2);
        book(_p["anchg"],                1, 1, 3);
        book(_p["pmaxnchg"],             2, 1, 1);
        book(_p["pminnchg"],             2, 1, 2);
        book(_p["pdifnchg"],             2, 1, 3);
        book(_p["tcptsum"],              3, 1, 1);
        book(_p["pcptsum"],              3, 1, 2);
        book(_p["acptsum"],              3, 1, 3);
        book(_p["pmaxcptsum"],           4, 1, 1);
        book(_p["pmincptsum"],           4, 1, 2);
        book(_p["pdifcptsum"],           4, 1, 3);
        book(_p["tcptave"],              5, 1, 1);
        book(_p["pcptave"],              5, 1, 2);
        book(_p["tcptmax"],              6, 1, 1);
        book(_p["pcptmax"],              6, 1, 2);
        book(_p["zptvsnchg"],            7, 1, 1);
        book(_p["cptavevsnchg"],         8, 1, 1);
        book(_p["cptavevsnchgsmallzpt"], 9, 1, 1);
      }

      if (_mode == 0 || _mode == 2) {
        book(_p["tnchg"],      10, 1, 1);
        book(_p["pnchg"],      10, 1, 2);
        book(_p["anchg"],      10, 1, 3);
        book(_p["pmaxnchg"],   11, 1, 1);
        book(_p["pminnchg"],   11, 1, 2);
        book(_p["pdifnchg"],   11, 1, 3);
        book(_p["tcptsum"],    12, 1, 1);
        book(_p["pcptsum"],    12, 1, 2);
        book(_p["acptsum"],    12, 1, 3);
        book(_p["pmaxcptsum"], 13, 1, 1);
        book(_p["pmincptsum"], 13, 1, 2);
        book(_p["pdifcptsum"], 13, 1, 3);
        book(_p["pcptave"],    14, 1, 1);
        book(_p["pcptmax"],    15, 1, 1);
      }
    }


    /// Do the analysis
    void analyze(const Event& event) {

      if (_mode == 0 || _mode == 1)  doDYanalysis(event);
      if (_mode == 0 || _mode == 2)  doQCDanalysis(event);

    }


    void doDYanalysis(const Event& e) {

      const FinalState& fs = apply<FinalState>(e, "CFS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 1) {
        MSG_DEBUG("Failed multiplicity cut");
        vetoEvent;
      }

      // Get the leptons
      const Particles& leptons = apply<ChargedLeptons>(e, "CL").chargedLeptons();

      // We want exactly two leptons of the same flavour.
      MSG_DEBUG("lepton multiplicity = " << leptons.size());
      if (leptons.size() != 2 || leptons[0].pid() != -leptons[1].pid() ) vetoEvent;

      // Lepton pT > 20 GeV
      if (leptons[0].pT()/GeV <= 20 || leptons[1].pT()/GeV <= 20) vetoEvent;

      // Lepton pair should have an invariant mass between 70 and 110 and |eta| < 6
      const FourMomentum dilepton = leptons[0].momentum() + leptons[1].momentum();
      if (!inRange(dilepton.mass()/GeV, 70., 110.) || fabs(dilepton.eta()) >= 6) vetoEvent;
      MSG_DEBUG("Dilepton mass = " << dilepton.mass()/GeV << " GeV");
      MSG_DEBUG("Dilepton pT   = " << dilepton.pT()/GeV << " GeV");

      // Calculate the observables
      size_t   numToward(0),     numAway(0);
      long int numTrans1(0),     numTrans2(0);
      double ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
      double ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);
      const double phiZ = dilepton.azimuthalAngle();
      const double pTZ  = dilepton.pT();
      /// @todo Replace with for
      for (Particles::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
        // Don't use the leptons
        /// @todo Replace with PID::isLepton
        if (abs(p->pid()) < 20) continue;

        const double dPhi = deltaPhi(p->momentum().phi(), phiZ);
        const double pT = p->pT();
        double rotatedphi = p->momentum().phi() - phiZ;
        while (rotatedphi < 0) rotatedphi += 2*PI;

        if (dPhi < PI/3.0) {
          ptSumToward += pT;
          ++numToward;
          if (pT > ptMaxToward)
            ptMaxToward = pT;
        } else if (dPhi < 2*PI/3.0) {
          if (rotatedphi <= PI) {
            ptSumTrans1 += pT;
            ++numTrans1;
            if (pT > ptMaxTrans1)
              ptMaxTrans1 = pT;
          }
          else {
            ptSumTrans2 += pT;
            ++numTrans2;
            if (pT > ptMaxTrans2)
              ptMaxTrans2 = pT;
          }
        } else {
          ptSumAway += pT;
          ++numAway;
          if (pT > ptMaxAway)
            ptMaxAway = pT;
        }
        // We need to subtract the two leptons from the number of particles to get the correct multiplicity
        _p["cptavevsnchg"]->fill(numParticles-2, pT);
        if (pTZ < 10)
          _p["cptavevsnchgsmallzpt"]->fill(numParticles-2, pT);
      }

      // Fill the histograms
      _p["tnchg"]->fill(pTZ, numToward/(4*PI/3));
      _p["pnchg"]->fill(pTZ, (numTrans1+numTrans2)/(4*PI/3));
      _p["pmaxnchg"]->fill(pTZ, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _p["pminnchg"]->fill(pTZ, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _p["pdifnchg"]->fill(pTZ, abs(numTrans1-numTrans2)/(2*PI/3));
      _p["anchg"]->fill(pTZ, numAway/(4*PI/3));

      _p["tcptsum"]->fill(pTZ, ptSumToward/(4*PI/3));
      _p["pcptsum"]->fill(pTZ, (ptSumTrans1+ptSumTrans2)/(4*PI/3));
      _p["pmaxcptsum"]->fill(pTZ, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3));
      _p["pmincptsum"]->fill(pTZ, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3));
      _p["pdifcptsum"]->fill(pTZ, fabs(ptSumTrans1-ptSumTrans2)/(2*PI/3));
      _p["acptsum"]->fill(pTZ, ptSumAway/(4*PI/3));

      if (numToward > 0) {
        _p["tcptave"]->fill(pTZ, ptSumToward/numToward);
        _p["tcptmax"]->fill(pTZ, ptMaxToward);
      }
      if ((numTrans1+numTrans2) > 0) {
        _p["pcptave"]->fill(pTZ, (ptSumTrans1+ptSumTrans2)/(numTrans1+numTrans2));
        _p["pcptmax"]->fill(pTZ, (ptMaxTrans1 > ptMaxTrans2 ? ptMaxTrans1 : ptMaxTrans2));
      }

      // We need to subtract the two leptons from the number of particles to get the correct multiplicity
      _p["zptvsnchg"]->fill(numParticles-2, pTZ);
    }


    void doQCDanalysis(const Event& e) {

      const FinalState& fsj = apply<FinalState>(e, "FSJ");
      if (fsj.particles().size() < 1) {
        MSG_DEBUG("Failed multiplicity cut");
        vetoEvent;
      }

      const Jets& jets = apply<FastJets>(e, "MidpointJets").jetsByPt();
      MSG_DEBUG("Jet multiplicity = " << jets.size());

      // We require the leading jet to be within |eta|<2
      if (jets.size() < 1 || fabs(jets[0].eta()) >= 2) {
        MSG_DEBUG("Failed leading jet cut");
        vetoEvent;
      }

      const double jetphi = jets[0].phi();
      const double jeteta = jets[0].eta();
      const double jetpT  = jets[0].pT();
      MSG_DEBUG("Leading jet: pT = " << jetpT
                << ", eta = " << jeteta << ", phi = " << jetphi);

      // Get the final states to work with for filling the distributions
      const FinalState& cfs = apply<ChargedFinalState>(e, "CFS");

      size_t numToward(0), numAway(0);
      long int numTrans1(0),     numTrans2(0);
      double ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
      double ptMaxOverall(0.0), ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);

      // Calculate all the charged stuff
      for (const Particle& p : cfs.particles()) {
        const double dPhi = deltaPhi(p.phi(), jetphi);
        const double pT = p.pT();
        const double phi = p.phi();
        double rotatedphi = phi - jetphi;
        while (rotatedphi < 0) rotatedphi += 2*PI;

        if (pT > ptMaxOverall) {
          ptMaxOverall = pT;
        }

        if (dPhi < PI/3.0) {
          ptSumToward += pT;
          ++numToward;
          if (pT > ptMaxToward) ptMaxToward = pT;
        }
        else if (dPhi < 2*PI/3.0) {
          if (rotatedphi <= PI) {
            ptSumTrans1 += pT;
            ++numTrans1;
            if (pT > ptMaxTrans1) ptMaxTrans1 = pT;
          } else {
            ptSumTrans2 += pT;
            ++numTrans2;
            if (pT > ptMaxTrans2) ptMaxTrans2 = pT;
          }
        }
        else {
          ptSumAway += pT;
          ++numAway;
          if (pT > ptMaxAway) ptMaxAway = pT;
        }
      } // end charged particle loop

      // Fill the histograms
      _p["tnchg"]->fill(jetpT/GeV, numToward/(4*PI/3));
      _p["pnchg"]->fill(jetpT/GeV, (numTrans1+numTrans2)/(4*PI/3));
      _p["pmaxnchg"]->fill(jetpT/GeV, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _p["pminnchg"]->fill(jetpT/GeV, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _p["pdifnchg"]->fill(jetpT/GeV, abs(numTrans1-numTrans2)/(2*PI/3));
      _p["anchg"]->fill(jetpT/GeV, numAway/(4*PI/3));

      _p["tcptsum"]->fill(jetpT/GeV, ptSumToward/GeV/(4*PI/3));
      _p["pcptsum"]->fill(jetpT/GeV, (ptSumTrans1+ptSumTrans2)/GeV/(4*PI/3));
      _p["pmaxcptsum"]->fill(jetpT/GeV, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/GeV/(2*PI/3));
      _p["pmincptsum"]->fill(jetpT/GeV, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/GeV/(2*PI/3));
      _p["pdifcptsum"]->fill(jetpT/GeV, fabs(ptSumTrans1-ptSumTrans2)/GeV/(2*PI/3));
      _p["acptsum"]->fill(jetpT/GeV, ptSumAway/GeV/(4*PI/3));

      if ((numTrans1+numTrans2) > 0) {
        _p["pcptave"]->fill(jetpT/GeV, (ptSumTrans1+ptSumTrans2)/GeV/(numTrans1+numTrans2));
        _p["pcptmax"]->fill(jetpT/GeV, (ptMaxTrans1 > ptMaxTrans2 ? ptMaxTrans1 : ptMaxTrans2)/GeV);
      }
    }


    // void finalize() {    }

    /// @}


  private:

    size_t _mode;
    map<string, Profile1DPtr> _p;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2010_I849042, CDF_2010_S8591881);

}
