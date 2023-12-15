// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/GammaGammaKinematics.hh"

namespace Rivet {


  /// @brief OPAL charged hadron production in photon-photon collisions.
  class OPAL_2007_I734955 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2007_I734955);

    /// @name Analysis methods
    /// @{

    // Book projections and histograms
    void init() {

      // Kinematics for leptons and intermediate photons.
      const GammaGammaKinematics gmgmkin;
      declare(gmgmkin, "Kinematics");

      const GammaGammaFinalState gmgmfs(Cuts::abseta < 1.5 && Cuts::pT > 0.12*GeV && Cuts::charge != 0, gmgmkin);
      declare(gmgmfs, "FS");

      // pT spectra of charged particles in W bins.
      book(_h_pTch[0], 1, 1, 1);
      book(_h_pTch[1], 2, 1, 1);
      book(_h_pTch[2], 3, 1, 1);
      book(_h_pTch[3], 4, 1, 1);

      // eta distributions of charged particles in W bins.
      book(_h_etach[0], 5, 1, 1);
      book(_h_etach[1], 6, 1, 1);
      book(_h_etach[2], 7, 1, 1);
      book(_h_etach[3], 8, 1, 1);

      // pT spectra of charged particles with L3-type W cuts.
      book(_h_pTch2[0], 9, 1, 1);
      book(_h_pTch2[1], 10, 1, 1);

      // pT spectra of charged pions with L3-type W cuts.
      book(_h_pTpi[0], 11, 1, 1);
      book(_h_pTpi[1], 12, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& event) {

      // Calculate invariant mass of photon pair.
      const GammaGammaKinematics& kinematics = apply<GammaGammaKinematics>(event, "Kinematics");

      const double Wgmgm = sqrt(kinematics.W2());
      const double Q2A = kinematics.Q2().first;
      const double Q2B = kinematics.Q2().second;

      // Reject events where photon virtuality is too large.
      if ( Q2A > 4.5*GeV2 || Q2B > 4.5*GeV2 ) vetoEvent;

      // Loop through charged particles.
      const Particles& chParticles = apply<GammaGammaFinalState>(event, "FS").particles();
      for (const Particle& p : chParticles) {
        const double pT  = p.pT() / GeV;
        const double eta = p.abseta();

        // Fill charged-hadron pT-spectra.
        if ( inRange(Wgmgm, 10., 30. ) ) _h_pTch[0]->fill(pT);
        if ( inRange(Wgmgm, 30., 50. ) ) _h_pTch[1]->fill(pT);
        if ( inRange(Wgmgm, 50., 125.) ) _h_pTch[2]->fill(pT);
        if ( inRange(Wgmgm, 10., 125.) ) _h_pTch[3]->fill(pT);

        // Fill charged-hadron eta distributions.
        if (pT > 3.5) {
          if ( inRange(Wgmgm, 10., 30. ) ) _h_etach[0]->fill(eta);
          if ( inRange(Wgmgm, 30., 50. ) ) _h_etach[1]->fill(eta);
          if ( inRange(Wgmgm, 50., 125.) ) _h_etach[2]->fill(eta);
          if ( inRange(Wgmgm, 10., 125.) ) _h_etach[3]->fill(eta);
        }

        // Fill charged-hadron pT-spectra with only lower W bound.
        if (Wgmgm > 30.) _h_pTch2[0]->fill(pT);
        if (Wgmgm > 50.) _h_pTch2[1]->fill(pT);

        // Fill charged-pion pT-spectra with only lower W bound.
        if ((p.abspid() == PID::PIPLUS || p.abspid() == PID::PIMINUS) && eta < 1.0) {
          if (Wgmgm > 30.) _h_pTpi[0]->fill(pT);
          if (Wgmgm > 50.) _h_pTpi[1]->fill(pT);
        }

      }

    }


    // Finalize
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();
      scale(_h_pTch,  sf);
      scale(_h_etach, sf);
      scale(_h_pTch2, sf);
      scale(_h_pTpi,  sf);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_pTch[4], _h_etach[4], _h_pTch2[2], _h_pTpi[2];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(OPAL_2007_I734955);

}
