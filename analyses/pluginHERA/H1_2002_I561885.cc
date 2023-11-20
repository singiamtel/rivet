// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of D*+- meson production in deep inelastic scattering at HERA (H1)
  class H1_2002_I561885 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2002_I561885);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {


      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "Dstars");
      //Cuts::abspid == PID::DSTARPLUS

      // Initialise and register projections


      Histo1DPtr dummy; //Introducing


      // Book histograms
      book(_h["W_GeV"], 2, 1, 1);
      book(_h["p_tD*"], 3, 1, 1);
      book(_h["logx"], 4, 1, 1);
      book(_h["etaD*"], 5, 1, 1);
      book(_h["Q2"], 6, 1, 1);
      book(_h["Z_D*"], 7, 1, 1);

      book(_h_Q2eta, {-1.5, -0.5, 0.5, 1.5}, {"d08-x01-y01", "d08-x01-y02", "d08-x01-y03"});
      book(_h_Q2pt, {1.5, 4., 10.}, {"d09-x01-y01", "d09-x01-y02"});
      book(_h_eta1, {0., 0.25, 0.5, 1.}, {"d10-x01-y01", "d10-x01-y02", "d10-x01-y03"});
      book(_h_eta2, {1.5, 2.5, 4., 10.}, {"d11-x01-y01", "d11-x01-y02", "d11-x01-y03"});
      book(_h_zD1, {1.5, 2.5, 4., 10.}, {"d12-x01-y01", "d12-x01-y02", "d12-x01-y03"});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");

       // Q2 and inelasticity cuts
      if (!inRange(kin.Q2(), 1.0*GeV2, 100*GeV2)) vetoEvent;
      if (!inRange(kin.y(), 0.05, 0.7)) vetoEvent;

        // D* reconstruction
      // const Particles unstables = apply<ParticleFinder>(event, "Dstars").particles(Cuts::pT > 1.5*GeV && Cuts::abseta < 1.5);
      const Particles unstables = apply<ParticleFinder>(event, "Dstars").particles(Cuts::pT > 1.5*GeV );
      const Particles dstars = select(unstables, [](const Particle& p){ return p.abspid() == PID::DSTARPLUS; });
      if (dstars.empty()) vetoEvent;
      // MSG_DEBUG("#D* = " << dstars.size());
      //const Particle& dstar = dstars.front();
      for (const Particle& dstar : dstars){
        const double zD = (dstar.E() - dstar.pz()) / (2*kin.beamLepton().E()*kin.y());


        // Single-differential histograms
        _h["p_tD*"]->fill(dstar.pT()/GeV);
        _h["etaD*"]->fill(dstar.eta());
        _h["Z_D*"]->fill(zD/GeV);
        _h["Q2"]->fill(kin.Q2()/GeV2);
        _h["W_GeV"]->fill(sqrt(kin.W2())/GeV);
        _h["logx"]->fill(log10(kin.x()));

        // Double-differential (y,Q2) histograms
        _h_Q2eta->fill(dstar.eta(),kin.Q2());
        _h_Q2pt->fill(dstar.pT(),kin.Q2());
        _h_eta1->fill(zD, dstar.eta());
        _h_eta2->fill(dstar.pT(), dstar.eta());
        _h_zD1->fill(dstar.pT(), zD/GeV);

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


     const double sf = crossSection()/nanobarn/sumOfWeights();
      scale(_h["p_tD*"], sf);
      scale(_h["etaD*"], sf);
      scale(_h["Z_D*"], sf);
      scale(_h["Q2"], sf);
      scale(_h["W_GeV"], sf);
      scale(_h["logx"], sf);

      scale(_h_Q2eta, sf);
      scale(_h_Q2pt, sf);
      scale(_h_eta1, sf);
      scale(_h_eta2, sf);
      scale(_h_zD1, sf);


    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    Histo1DGroupPtr _h_eta1,_h_eta2,_h_zD1, _h_Q2eta, _h_Q2pt;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_2002_I561885);

}
