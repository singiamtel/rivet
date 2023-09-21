#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Drell-Yan dimuon absolute cross-sections in 800 GeV pp and pd collisions
  class NUSEA_2003_I613362 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NUSEA_2003_I613362);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState fs;
      declare(fs, "FS");
      Cut cut = Cuts::etaIn(-10.,10.);
      ZFinder zfinder(fs, cut, PID::MUON, 4.0*GeV, 100.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE );
      declare(zfinder, "ZFinder");

      // Booking histograms
      // hydrogen d01-d16
      book(_hist_M_xF, {-0.05, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8});
      for (auto& b : _hist_M_xF->bins()) {
        book(b, b.index(), 1, 1);
      }

      // deuterium d17-d32

      // hydrogen d40
      book(_hist_pT_M, {4.2, 5.2, 6.2, 7.2, 8.7, 10.85, 12.85});
      _hist_pT_M->maskBin(5); int idx = 0;
      for (auto& b : _hist_pT_M->bins()) {
        book(b, 40, 1, ++idx);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double sqrts_tol = 10.;
      if (!isCompatibleWithSqrtS(38.8*GeV, sqrts_tol)) {
        MSG_ERROR("Incorrect beam energy used: " << sqrtS()/GeV);
        throw Error("Unexpected sqrtS ! Only 38.8 GeV is supported");
      }

      // Muons
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() <= 0) vetoEvent;

      double Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
      double Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
      double Zpl   = zfinder.bosons()[0].momentum().pz()/GeV;
      double ZE    = zfinder.bosons()[0].momentum().E();
      double xF = 2.*Zpl/sqrtS();

      // Filling dimuon mass in bins of xF
      _hist_M_xF->fill(xF, Zmass/GeV, Zmass*sqr(Zmass));

      // Filling pT in bins of Zmass
      if ( xF > -0.05 && xF <= 0.15 ) {
        // Include here all factors which are run-dependent for later scaling
        if (Zpt > 0) _hist_pT_M->fill(Zmass,Zpt, 1./2./Zpt*2.*ZE/sqrtS());
      }

      MSG_DEBUG("Dimuon pT = "<< Zpt<<"   Dimuon E = ");
      MSG_DEBUG("DiMuon mass " << Zmass/GeV);
      MSG_DEBUG("DiMuon pT "<<Zpt);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // xf bin width = 0.2, x-section in picobarn
      const double scalefactor=crossSection()/picobarn/(sumOfWeights() * M_PI *0.2 );
      scale(_hist_pT_M, scalefactor);

      // x-section is quoted in nanobarn
      scale(_hist_M_xF, crossSection()/nanobarn/sumOfWeights());

      divByGroupWidth({_hist_pT_M, _hist_M_xF});
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DGroupPtr _hist_pT_M, _hist_M_xF;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(NUSEA_2003_I613362);

}
