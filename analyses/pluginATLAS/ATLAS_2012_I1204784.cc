// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// ATLAS Z phi* measurement
  class ATLAS_2012_I1204784 : public Analysis {
    public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2012_I1204784);


    public:

      /// Book histograms and initialise projections before the run
      void init() {
        Cut cuts = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
        DileptonFinder zfinder_dressed_el(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
        declare(zfinder_dressed_el, "DileptonFinder_dressed_el");
        DileptonFinder zfinder_bare_el(91.2*GeV, 0.0, cuts && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
        declare(zfinder_bare_el, "DileptonFinder_bare_el");
        DileptonFinder zfinder_dressed_mu(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::MUON, Cuts::massIn(66*GeV, 116*GeV));
        declare(zfinder_dressed_mu, "DileptonFinder_dressed_mu");
        DileptonFinder zfinder_bare_mu(91.2*GeV, 0.0, cuts && Cuts::abspid == PID::MUON, Cuts::massIn(66*GeV, 116*GeV));
        declare(zfinder_bare_mu, "DileptonFinder_bare_mu");

        // Book histograms
        // Single-differential plots
        book(_hist_zphistar_el_bare ,1, 1, 1);
        book(_hist_zphistar_mu_bare ,1, 1, 2);
        book(_hist_zphistar_el_dressed ,2, 1, 1);
        book(_hist_zphistar_mu_dressed ,2, 1, 2);

        // Double-differential plots
        book(_h_phistar_el_bare,    {0., 0.8, 1.6, 10.}, {"d03-x01-y01", "d03-x01-y02", "d03-x01-y03"});
        book(_h_phistar_el_dressed, {0., 0.8, 1.6, 10.}, {"d03-x02-y01", "d03-x02-y02", "d03-x02-y03"});
        book(_h_phistar_mu_bare,    {0., 0.8, 1.6, 10.}, {"d04-x01-y01", "d04-x01-y02", "d04-x01-y03"});
        book(_h_phistar_mu_dressed, {0., 0.8, 1.6, 10.}, {"d04-x02-y01", "d04-x02-y02", "d04-x02-y03"});

      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        const DileptonFinder& zfinder_dressed_el = apply<DileptonFinder>(event, "DileptonFinder_dressed_el");
        fillPlots(zfinder_dressed_el, _hist_zphistar_el_dressed, _h_phistar_el_dressed);

        const DileptonFinder& zfinder_bare_el = apply<DileptonFinder>(event, "DileptonFinder_bare_el");
        fillPlots(zfinder_bare_el, _hist_zphistar_el_bare, _h_phistar_el_bare);

        const DileptonFinder& zfinder_dressed_mu = apply<DileptonFinder>(event, "DileptonFinder_dressed_mu");
        fillPlots(zfinder_dressed_mu, _hist_zphistar_mu_dressed, _h_phistar_mu_dressed);

        const DileptonFinder& zfinder_bare_mu = apply<DileptonFinder>(event, "DileptonFinder_bare_mu");
        fillPlots(zfinder_bare_mu, _hist_zphistar_mu_bare, _h_phistar_mu_bare);
      }


      void fillPlots(const DileptonFinder& zfind, Histo1DPtr hist, Histo1DGroupPtr& binnedHist) {
        if (zfind.bosons().size() != 1) return;
        Particles leptons = sortBy(zfind.constituents(), cmpMomByPt);

        const FourMomentum lminus = leptons[0].charge() < 0 ? leptons[0].momentum() : leptons[1].momentum();
        const FourMomentum lplus = leptons[0].charge() < 0 ? leptons[1].momentum() : leptons[0].momentum();

        const double phi_acop = M_PI - deltaPhi(lminus, lplus);
        const double costhetastar = tanh((lminus.eta()-lplus.eta())/2.0);
        const double sin2thetastar = (costhetastar <= 1) ? 1.0 - sqr(costhetastar) : 0;
        const double phistar = tan(phi_acop/2.0) * sqrt(sin2thetastar);
        hist->fill(phistar);

        binnedHist->fill(zfind.bosons()[0].absrap(), phistar);
      }


      /// Normalise histograms etc., after the run
      void finalize() {
        normalize(_hist_zphistar_el_dressed);
        normalize(_hist_zphistar_el_bare);
        normalize(_hist_zphistar_mu_dressed);
        normalize(_hist_zphistar_mu_bare);

        normalize(_h_phistar_mu_dressed);
        normalize(_h_phistar_mu_bare);
        normalize(_h_phistar_el_bare);
        normalize(_h_phistar_el_dressed);
      }

      //@}


    private:

      Histo1DGroupPtr _h_phistar_mu_dressed;
      Histo1DGroupPtr _h_phistar_mu_bare;
      Histo1DGroupPtr _h_phistar_el_dressed;
      Histo1DGroupPtr _h_phistar_el_bare;

      Histo1DPtr _hist_zphistar_el_dressed;
      Histo1DPtr _hist_zphistar_el_bare;

      Histo1DPtr _hist_zphistar_mu_dressed;
      Histo1DPtr _hist_zphistar_mu_bare;

      Histo1DPtr _hist_zphistar_el_bare_1;
      Histo1DPtr _hist_zphistar_el_bare_2;
      Histo1DPtr _hist_zphistar_el_bare_3;

      Histo1DPtr _hist_zphistar_el_dressed_1;
      Histo1DPtr _hist_zphistar_el_dressed_2;
      Histo1DPtr _hist_zphistar_el_dressed_3;

      Histo1DPtr _hist_zphistar_mu_bare_1;
      Histo1DPtr _hist_zphistar_mu_bare_2;
      Histo1DPtr _hist_zphistar_mu_bare_3;

      Histo1DPtr _hist_zphistar_mu_dressed_1;
      Histo1DPtr _hist_zphistar_mu_dressed_2;
      Histo1DPtr _hist_zphistar_mu_dressed_3;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2012_I1204784);

}
