// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "fastjet/SISConePlugin.hh"

namespace Rivet {


  /// @brief Forward jet production in deep inelastic scattering at HERA (ZEUS)
  class ZEUS_1999_I470499 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1999_I470499);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs;
      declare(fs, "FS");

      // declare jets
      double jet_radius = 1.0;
      declare(FastJets(fs, JetAlg::PXCONE, jet_radius), "Jets");

      // declare DIS Kinematics
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");

      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["xbj"], 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      const FinalState& fs = apply<FinalState>(event, "FS");

      const size_t numParticles = fs.particles().size();

      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }

      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double xbj  = dk.x();
      double ybj = dk.y();
      double Q2 = dk.Q2()/GeV;

      //cut on y
      if (ybj<0.1) vetoEvent;
      if (4.5*pow(10,-4) > xbj && xbj > 4.5*pow(10,-2)) vetoEvent;

      //Frame transfer
      const LorentzTransform breitboost = dk.boostBreit();

      //on scattered lepton

      FourMomentum leptonMom = dl.out().mom();
      double enel = leptonMom.E();

      bool cut = enel>10*GeV;
      if (!cut) vetoEvent;

      //scattered jets
      const Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 5*GeV && Cuts::eta<2.6, cmpMomByEt);

      bool loopjet=false;
      for (const Jet& j : jets) {
        //cout << " j.pz "  << j.pz()/(820*GeV) << endl;
        if (j.pz()/(820*GeV)<0.036) continue;
        // cout << pow(j.Et(),2) << endl;

        if (0.5 > sqr(j.Et())/Q2 || sqr(j.Et())/Q2 > 2.0) continue;
        FourMomentum breMom = breitboost.transform(j.mom());
        //cout<< " pz " << breMom.pz() << endl;
        if (breMom.pz()<0) continue;

        loopjet=true;
      }
      // cout<< " loojet " << loopjet << endl;
      if (loopjet)  _h["xbj"]->fill(xbj);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h["xbj"], crossSection()/nanobarn/sumW());
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}


  };

  RIVET_DECLARE_PLUGIN(ZEUS_1999_I470499);
}
