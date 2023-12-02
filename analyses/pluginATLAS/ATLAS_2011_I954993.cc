// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief WZ fiducial cross-section measurement
  class ATLAS_2011_I954993 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I954993);


    /// @name Analysis methods
    /// @{

    /// Projection and histogram setup
    void init() {
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;
      DileptonFinder zfinder_e(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::ELECTRON,
                               Cuts::massIn(81.1876*GeV, 101.1876*GeV));
      declare(zfinder_e, "DileptonFinder_e");
      DileptonFinder zfinder_mu(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::MUON,
                                Cuts::massIn(81.1876*GeV, 101.1876*GeV));
      declare(zfinder_mu, "DileptonFinder_mu");

      declare("MET", MissingMomentum());
      VetoedFinalState weinput;
      weinput.addVetoOnThisFinalState(zfinder_e);
      LeptonFinder ef(weinput, 0.1, cuts && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");
      VetoedFinalState wminput;
      wminput.addVetoOnThisFinalState(zfinder_mu);
      LeptonFinder mf(wminput, 0.1, cuts && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      // Histograms
      book(_h_fiducial ,1,1,1);
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // Looking for a Z, exit if not found
      const DileptonFinder& zfinder_e = apply<DileptonFinder>(event, "DileptonFinder_e");
      const DileptonFinder& zfinder_mu = apply<DileptonFinder>(event, "DileptonFinder_mu");
      if (zfinder_e.bosons().size() != 1 && zfinder_mu.bosons().size() != 1) {
        MSG_DEBUG("No Z boson found, vetoing event");
        vetoEvent;
      }

      // Looking for a W, exit if not found
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.Et() < 25*GeV) vetoEvent;
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int iefound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);
      const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
      const int imfound = closestMatchIndex(mus, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);
      if (iefound < 0 && imfound < 0) {
        MSG_DEBUG("No W boson found, vetoing event");
        vetoEvent;
      }

      // If we find a W, make fiducial acceptance cuts and exit if not found
      if (iefound >= 0) {
        const Particle& e = es[iefound];
        const double mt = mT(pmiss, e);
        if (e.pT() < 20*GeV || mt < 20*GeV) vetoEvent;
      } else if (imfound >= 0) {
        const Particle& m = mus[imfound];
        const double mt = mT(pmiss, m);
        if (m.pT() < 20*GeV || mt < 20*GeV) vetoEvent;
      }

      // Update the fiducial cross-section histogram
      /// @todo Would be better as a counter
      _h_fiducial->fill(7000);
    }


    /// Finalize
    void finalize() {
      scale(_h_fiducial, crossSection()/femtobarn/sumOfWeights());
    }

    /// @}


  private:

    /// Histogram
    Histo1DPtr _h_fiducial;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2011_I954993);

}
