// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief CMS Z pT and rapidity in Drell-Yan events at 7 TeV
  ///
  /// @author Justin Hugon, Luca Perrozzi
  class CMS_2012_I941555 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I941555);


    /// @name Analysis methods
    /// @{

    void init() {

      // Set up projections
      /// @todo Really?: was DileptonFinder zfinder_dressed_mu_pt(-2.1, 2.1, 20, PID::MUON, 60*GeV, 120*GeV, 0.2, false, true);
      Cut cuts = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;
      DileptonFinder zfinder_dressed_mu_pt(91.2*GeV, 0.2, cuts && Cuts::abspid == PID::MUON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zfinder_dressed_mu_pt, "DileptonFinder_dressed_mu_pt");
      DileptonFinder zfinder_dressed_el_pt(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::ELECTRON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zfinder_dressed_el_pt, "DileptonFinder_dressed_el_pt");

      DileptonFinder zfinder_dressed_mu_rap(91.2*GeV, 0.1, Cuts::abspid == PID::MUON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zfinder_dressed_mu_rap, "DileptonFinder_dressed_mu_rap");
      DileptonFinder zfinder_dressed_el_rap(91.2*GeV, 0.1, Cuts::abspid == PID::ELECTRON, Cuts::massIn(60*GeV, 120*GeV));
      declare(zfinder_dressed_el_rap, "DileptonFinder_dressed_el_rap");

      // Book histograms
      book(_hist_zrap_mu_dressed      , 1, 1, 1);  // muon "dressed" rapidity
      book(_hist_zrap_el_dressed      , 1, 1, 2);  // electron "dressed" rapidity
      book(_hist_zrap_comb_dressed    , 1, 1, 3);  // electron "dressed" rapidity

      book(_hist_zpt_mu_dressed       , 2, 1, 1);  // muon "dressed" pt
      book(_hist_zpt_el_dressed       , 2, 1, 2);  // electron "dressed" pt
      book(_hist_zpt_comb_dressed     , 2, 1, 3);  // electron "dressed" pt

      book(_hist_zptpeak_mu_dressed   , 3, 1, 1);  // muon "dressed" pt peak
      book(_hist_zptpeak_el_dressed   , 3, 1, 2);  // electron "dressed" pt peak
      book(_hist_zptpeak_comb_dressed , 3, 1, 3);  // electron "dressed" pt peak
    }


    /// Do the analysis
    void analyze(const Event& evt) {

      const DileptonFinder& zfinder_dressed_mu_rap = apply<DileptonFinder>(evt, "DileptonFinder_dressed_mu_rap");
      if (!zfinder_dressed_mu_rap.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_mu_rap.bosons()[0].momentum();
        _hist_zrap_mu_dressed->fill(pZ.absrapidity());
        _hist_zrap_comb_dressed->fill(pZ.absrapidity());
      }

      const DileptonFinder& zfinder_dressed_mu_pt = apply<DileptonFinder>(evt, "DileptonFinder_dressed_mu_pt");
      if (!zfinder_dressed_mu_pt.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_mu_pt.bosons()[0].momentum();
        _hist_zpt_mu_dressed->fill(pZ.pT()/GeV);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV);
        if (pZ.pT() < 30*GeV) {
          _hist_zptpeak_mu_dressed->fill(pZ.pT()/GeV);
          _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV);
        }
      }

      const DileptonFinder& zfinder_dressed_el_rap = apply<DileptonFinder>(evt, "DileptonFinder_dressed_el_rap");
      if (!zfinder_dressed_el_rap.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_el_rap.bosons()[0].momentum();
        _hist_zrap_el_dressed->fill(abs(pZ.rapidity()/GeV));
        _hist_zrap_comb_dressed->fill(abs(pZ.rapidity()/GeV));
      }

      const DileptonFinder& zfinder_dressed_el_pt = apply<DileptonFinder>(evt, "DileptonFinder_dressed_el_pt");
      if (!zfinder_dressed_el_pt.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_el_pt.bosons()[0].momentum();
        _hist_zpt_el_dressed->fill(pZ.pT()/GeV);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV);
        if (pZ.pT() < 30*GeV) {
          _hist_zptpeak_el_dressed->fill(pZ.pT()/GeV);
          _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV);
        }
      }

    }


    void finalize() {
      normalize(_hist_zrap_mu_dressed);
      normalize(_hist_zpt_mu_dressed);
      normalize(_hist_zptpeak_mu_dressed);

      normalize(_hist_zrap_el_dressed);
      normalize(_hist_zpt_el_dressed);
      normalize(_hist_zptpeak_el_dressed);

      normalize(_hist_zrap_comb_dressed);
      normalize(_hist_zpt_comb_dressed);
      normalize(_hist_zptpeak_comb_dressed);
    }

    /// @}


  private:

    Histo1DPtr _hist_zrap_mu_dressed;
    Histo1DPtr _hist_zpt_mu_dressed;
    Histo1DPtr _hist_zptpeak_mu_dressed;

    Histo1DPtr _hist_zrap_el_dressed;
    Histo1DPtr _hist_zpt_el_dressed;
    Histo1DPtr _hist_zptpeak_el_dressed;

    Histo1DPtr _hist_zrap_comb_dressed;
    Histo1DPtr _hist_zpt_comb_dressed;
    Histo1DPtr _hist_zptpeak_comb_dressed;

  };


  RIVET_DECLARE_PLUGIN(CMS_2012_I941555);

}
