// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// Precise study of Z pT using novel phi* technique
  class D0_2010_S8821313 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_S8821313);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      FinalState fs;
      Cut cuts_e = (Cuts::abseta < 1.1 || Cuts::absetaIn( 1.5,  3.0)) && Cuts::pT > 20*GeV;
      DileptonFinder zfinder_ee(91.2*GeV, 0.2, cuts_e && Cuts::abspid == PID::ELECTRON, Cuts::massIn(70*GeV, 110*GeV));
      declare(zfinder_ee, "zfinder_ee");
      Cut cuts_m = Cuts::abseta < 2 && Cuts::pT > 15*GeV;
      DileptonFinder zfinder_mm(91.2*GeV, 0.0, cuts_m && Cuts::abspid == PID::MUON, Cuts::massIn(70*GeV, 110*GeV));
      declare(zfinder_mm, "zfinder_mm");

      /// Book histograms here
      book(_h_phistar_ee, {0., 1., 2., 10.}, {"d01-x01-y01", "d01-x01-y02", "d01-x01-y03"});
      book(_h_phistar_mm, {0., 1., 2.}, {"d02-x01-y01", "d02-x01-y02"});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DileptonFinder& zfinder_ee = apply<DileptonFinder>(event, "zfinder_ee");
      if (zfinder_ee.bosons().size() == 1) {
        Particles ee = zfinder_ee.constituents();
        std::sort(ee.begin(), ee.end(), cmpMomByPt);
        const FourMomentum& eminus = PID::charge3(ee[0].pid()) < 0 ? ee[0].momentum() : ee[1].momentum();
        const FourMomentum& eplus  = PID::charge3(ee[0].pid()) < 0 ? ee[1].momentum() : ee[0].momentum();
        double phi_acop = M_PI - mapAngle0ToPi(eminus.phi() - eplus.phi());
        double costhetastar = tanh((eminus.eta() - eplus.eta())/2);
        double sin2thetastar = 1 - sqr(costhetastar);
        if (sin2thetastar < 0) sin2thetastar = 0;
        const double phistar = tan(phi_acop/2) * sqrt(sin2thetastar);
        const FourMomentum& zmom = zfinder_ee.bosons()[0].momentum();
        _h_phistar_ee->fill(zmom.rapidity(), phistar);
      }

      const DileptonFinder& zfinder_mm = apply<DileptonFinder>(event, "zfinder_mm");
      if (zfinder_mm.bosons().size() == 1) {
        Particles mm = zfinder_mm.constituents();
        std::sort(mm.begin(), mm.end(), cmpMomByPt);
        const FourMomentum& mminus = PID::charge3(mm[0].pid()) < 0 ? mm[0].momentum() : mm[1].momentum();
        const FourMomentum& mplus  = PID::charge3(mm[0].pid()) < 0 ? mm[1].momentum() : mm[0].momentum();
        double phi_acop = M_PI - mapAngle0ToPi(mminus.phi() - mplus.phi());
        double costhetastar = tanh((mminus.eta() - mplus.eta())/2);
        double sin2thetastar = 1 - sqr(costhetastar);
        if (sin2thetastar < 0) sin2thetastar = 0;
        const double phistar = tan(phi_acop/2) * sqrt(sin2thetastar);
        const FourMomentum& zmom = zfinder_mm.bosons()[0].momentum();
        _h_phistar_mm->fill(zmom.rapidity(), phistar);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize({_h_phistar_ee, _h_phistar_mm});
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_phistar_ee, _h_phistar_mm;
    /// @}


  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_S8821313, D0_2010_I871787);

}
