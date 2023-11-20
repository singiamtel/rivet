// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  /// @brief Pseudorapidity distributions of charged particles at sqrt{s} = 0.9 and 7 TeV
  class CMS_2010_PAS_QCD_10_024 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2010_PAS_QCD_10_024);


    void init() {
      declare(ChargedFinalState((Cuts::etaIn(-0.8, 0.8) && Cuts::pT >= 0.5*GeV)), "CFS_08_05");
      declare(ChargedFinalState((Cuts::etaIn(-0.8, 0.8) && Cuts::pT >= 1.0*GeV)), "CFS_08_10");
      declare(ChargedFinalState((Cuts::etaIn(-2.4, 2.4) && Cuts::pT >= 0.5*GeV)), "CFS_24_05");
      declare(ChargedFinalState((Cuts::etaIn(-2.4, 2.4) && Cuts::pT >= 1.0*GeV)), "CFS_24_10");

      size_t offset = 0;
      if (isCompatibleWithSqrtS(7000*GeV)) offset = 0;
      if (isCompatibleWithSqrtS(900*GeV)) offset = 4;
      book(_hist_dNch_deta_pt05_eta08 ,1+offset, 1, 1);
      book(_hist_dNch_deta_pt10_eta08 ,2+offset, 1, 1);
      book(_hist_dNch_deta_pt05_eta24 ,3+offset, 1, 1);
      book(_hist_dNch_deta_pt10_eta24 ,4+offset, 1, 1);

      book(_weight_pt05_eta08, "_pt05_eta08");
      book(_weight_pt10_eta08, "_pt10_eta08");
      book(_weight_pt05_eta24, "_pt05_eta24");
      book(_weight_pt10_eta24, "_pt10_eta24");
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs_08_05 = apply<ChargedFinalState>(event, "CFS_08_05");
      const ChargedFinalState& cfs_08_10 = apply<ChargedFinalState>(event, "CFS_08_10");
      const ChargedFinalState& cfs_24_05 = apply<ChargedFinalState>(event, "CFS_24_05");
      const ChargedFinalState& cfs_24_10 = apply<ChargedFinalState>(event, "CFS_24_10");

      // Plot distributions
      if(!cfs_08_05.particles().empty()) _weight_pt05_eta08->fill();
      if(!cfs_24_05.particles().empty()) _weight_pt05_eta24->fill();
      for (const Particle& p : cfs_24_05.particles()) {
        _hist_dNch_deta_pt05_eta24->fill(p.eta());
        if(!cfs_08_05.particles().empty())
          _hist_dNch_deta_pt05_eta08->fill(p.eta());
      }
      if(!cfs_08_10.particles().empty()) _weight_pt10_eta08->fill();
      if(!cfs_24_10.particles().empty()) _weight_pt10_eta24->fill();
      for (const Particle& p : cfs_24_10.particles()) {
        _hist_dNch_deta_pt10_eta24->fill(p.eta());
        if(!cfs_08_10.particles().empty())
          _hist_dNch_deta_pt10_eta08->fill(p.eta());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_dNch_deta_pt05_eta08,1./ *_weight_pt05_eta08);
      scale(_hist_dNch_deta_pt10_eta08,1./ *_weight_pt10_eta08);
      scale(_hist_dNch_deta_pt05_eta24,1./ *_weight_pt05_eta24);
      scale(_hist_dNch_deta_pt10_eta24,1./ *_weight_pt10_eta24);
    }


  private:

    Histo1DPtr _hist_dNch_deta_pt05_eta08;
    Histo1DPtr _hist_dNch_deta_pt10_eta08;
    Histo1DPtr _hist_dNch_deta_pt05_eta24;
    Histo1DPtr _hist_dNch_deta_pt10_eta24;
    CounterPtr _weight_pt05_eta08,_weight_pt10_eta08,_weight_pt05_eta24,_weight_pt10_eta24;
  };


  RIVET_DECLARE_PLUGIN(CMS_2010_PAS_QCD_10_024);

}
