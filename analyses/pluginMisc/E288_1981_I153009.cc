// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Measurement of dimuon continuum in proton-nucleus collisions
  class E288_1981_I153009 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E288_1981_I153009);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;

      ZFinder zfinder(fs, Cuts::abseta < 15.0, PID::MUON, 3.5*GeV, 30.0*GeV, 0.1, PhotonOrigin::NONE );
      declare(zfinder, "ZFinder");

      // Book histograms
      // 400 GeV and y = 0.03
      book(_hist_pT_M_400, {5., 6., 7., 8., 9., 10., 11., 12., 13., 14.});
      for (auto& b : _hist_pT_M_400->bins()) {
        book(b, 9, 1, b.index());
      }

      int Nbin = 50;
      book(_h_m_DiMuon,  "DiMuon_mass", Nbin, 0.0,30.0);
      book(_h_pT_DiMuon, "DiMuon_pT",   Nbin, 0.0,20.0);
      book(_h_y_DiMuon,  "DiMuon_y",    Nbin,-8.0, 8.0);
      book(_h_xF_DiMuon, "DiMuon_xF",   Nbin,-1.5, 1.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double sqrts_tol = 10. ;
      if (!isCompatibleWithSqrtS(27.4*GeV, sqrts_tol)) {
        MSG_ERROR("Incorrect beam energy used: " << sqrtS()/GeV);
        throw Error("Unexpected sqrtS ! Only 27.4 GeV is supported");
      }

      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() >= 1) {

        double Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
        double Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
        double Zpl   = zfinder.bosons()[0].momentum().pz()/GeV;
        double Zy    = zfinder.bosons()[0].momentum().rapidity();
        //double ZE    = zfinder.bosons()[0].momentum().E();

        double xf = 2.*Zpl/sqrtS() ;
        _h_xF_DiMuon->fill(xf);
        _h_m_DiMuon->fill(Zmass/GeV);
        _h_pT_DiMuon->fill(Zpt);
        _h_y_DiMuon ->fill(Zy);
        double Zymin = -1.0;
        double Zymax = 1.03;
        double Z_y_width = Zymax - Zymin ;
        if ( Zy > Zymin && Zy < Zymax ) {
          // Edsigma^3/dp^3 = 2E/(pi*sqrts)dsigma/dx_F/dq_T^2 = 1/pi dsigma/dy/dq_T^2
          // normalisation of Zy bin width = Zwidth
          if (Zpt > 0) _hist_pT_M_400->fill(Zmass,Zpt, 1./2./Zpt/Z_y_width);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_pT_M_400, crossSection()/femtobarn/(sumOfWeights() * M_PI));
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _hist_pT_M_400;
    Histo1DPtr _h_m_DiMuon ;
    Histo1DPtr _h_pT_DiMuon;
    Histo1DPtr _h_y_DiMuon;
    Histo1DPtr _h_xF_DiMuon;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(E288_1981_I153009);

}
