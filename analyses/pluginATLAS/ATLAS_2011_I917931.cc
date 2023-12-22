// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief ATLAS Z pT in Drell-Yan events at 7 TeV
  ///
  /// @author Elena Yatsenko, Judith Katzy
  class ATLAS_2011_I917931 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I917931);


    /// @name Analysis methods
    /// @{

    void init() {

      // Set up projections
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      DileptonFinder zfinder_dressed_el(91.2*GeV, 0.1, cut && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_dressed_el, "DileptonFinder_dressed_el");
      DileptonFinder zfinder_bare_el(91.2*GeV, 0.0, cut && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_bare_el, "DileptonFinder_bare_el");
      DileptonFinder zfinder_dressed_mu(91.2*GeV, 0.1, cut && Cuts::abspid == PID::MUON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_dressed_mu, "DileptonFinder_dressed_mu");
      DileptonFinder zfinder_bare_mu(91.2*GeV, 0.0, cut && Cuts::abspid == PID::MUON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_bare_mu, "DileptonFinder_bare_mu");

      // Book histograms
      book(_hist_zpt_el_dressed     ,1, 1, 2);  // electron "dressed"
      book(_hist_zpt_el_bare        ,1, 1, 3);  // electron "bare"
      book(_hist_zpt_mu_dressed     ,2, 1, 2);  // muon "dressed"
      book(_hist_zpt_mu_bare        ,2, 1, 3);  // muon "bare"

      book(_sumw_el_bare, "_sumw_el_bare");
      book(_sumw_el_dressed, "_sumw_el_dressed");
      book(_sumw_mu_bare, "_sumw_mu_bare");
      book(_sumw_mu_dressed, "_sumw_mu_dressed");
    }


    /// Do the analysis
    void analyze(const Event& evt) {
      const DileptonFinder& zfinder_dressed_el = apply<DileptonFinder>(evt, "DileptonFinder_dressed_el");
      if (!zfinder_dressed_el.bosons().empty()) {
        _sumw_el_dressed->fill();
        const FourMomentum pZ = zfinder_dressed_el.bosons()[0].momentum();
        _hist_zpt_el_dressed->fill(pZ.pT()/GeV);
      }

      const DileptonFinder& zfinder_bare_el = apply<DileptonFinder>(evt, "DileptonFinder_bare_el");
      if (!zfinder_bare_el.bosons().empty()) {
        _sumw_el_bare->fill();
	    const FourMomentum pZ = zfinder_bare_el.bosons()[0].momentum();
        _hist_zpt_el_bare->fill(pZ.pT()/GeV);
      }

      const DileptonFinder& zfinder_dressed_mu = apply<DileptonFinder>(evt, "DileptonFinder_dressed_mu");
      if (!zfinder_dressed_mu.bosons().empty()) {
        _sumw_mu_dressed->fill();
        const FourMomentum pZ = zfinder_dressed_mu.bosons()[0].momentum();
        _hist_zpt_mu_dressed->fill(pZ.pT()/GeV);
      }

      const DileptonFinder& zfinder_bare_mu = apply<DileptonFinder>(evt, "DileptonFinder_bare_mu");
      if (!zfinder_bare_mu.bosons().empty()) {
        _sumw_mu_bare->fill();
        const FourMomentum pZ = zfinder_bare_mu.bosons()[0].momentum();
        _hist_zpt_mu_bare->fill(pZ.pT()/GeV);
      }

    }


    void finalize() {
      if (_sumw_el_dressed->val() != 0) scale(_hist_zpt_el_dressed, 1/ *_sumw_el_dressed);
      if (_sumw_el_bare->val()    != 0) scale(_hist_zpt_el_bare,    1/ *_sumw_el_bare);
      if (_sumw_mu_dressed->val() != 0) scale(_hist_zpt_mu_dressed, 1/ *_sumw_mu_dressed);
      if (_sumw_mu_bare->val()    != 0) scale(_hist_zpt_mu_bare,    1/ *_sumw_mu_bare);
    }

    /// @}


    private:

	CounterPtr _sumw_el_bare, _sumw_el_dressed;
	CounterPtr _sumw_mu_bare, _sumw_mu_dressed;

	Histo1DPtr _hist_zpt_el_dressed;
	Histo1DPtr _hist_zpt_el_bare;
	Histo1DPtr _hist_zpt_mu_dressed;
	Histo1DPtr _hist_zpt_mu_bare;
  };



  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_I917931, ATLAS_2011_S9131140);

}
