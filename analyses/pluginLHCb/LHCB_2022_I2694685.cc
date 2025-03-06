// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet
{

  /// @brief Prompt D0 production in pPb collisions at 8.16 TeV c.o.m. energy per nucleon
  class LHCB_2022_I2694685 : public Analysis
  {

  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2022_I2694685);

    void init()
    {
      // This booking may be optimised but since BinnedHistogram changes in Rivet 4
      // we keep it this way for pedagogical reasons.
      book(_pt_hist, {1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00,
        -2.75, -2.50, -3.00, -3.25, -3.50, -3.75, -4.00, -4.25, -4.50, -4.75, -5.00}, 
        {"d01-x01-y01", "d01-x01-y02", "d01-x01-y03", "d01-x01-y04", "d01-x01-y05", "d01-x01-y06", 
          "d01-x01-y07", "d01-x01-y08", "d01-x01-y09", "d01-x01-y10",
        "d02-x01-y01", "d02-x01-y02", "d02-x01-y03", "d02-x01-y04", "d02-x01-y05", "d02-x01-y06",
          "d02-x01-y07", "d02-x01-y08", "d02-x01-y09", "d02-x01-y10", "d02-x01-y10"});
      book(_forw_to_backw_rat, 5, 1, 1);

      // Define histograms with transverse-momentum intervals to compute
      // forward-to-backward ratio
      book(_forw_y_pt_hist, "TMP/_forw_y_pt_hist", refData(5, 1, 1));
      book(_backw_y_pt_hist, "TMP/_backw_y_pt_hist", refData(5, 1, 1));

      // Select D0 and D0bar mesons
      declare(UnstableParticles(
        (Cuts::abspid == PID::D0) &&
        (Cuts::pT < 30 * GeV) &&
        (Cuts::rapIn(1.5, 4.0) || Cuts::rapIn(-5.0, -2.5))), "UPs");
    }

    void analyze(const Event &ev)
    {
      // Apply UnstableParticles projection to get entire decay chain of every particle
      const UnstableParticles &unst_parts = apply<UnstableParticles>(ev, "UPs");
      for (const Particle &part : unst_parts.particles())
      {
        // veto particle if beauty parent
        if (part.fromBottom()) continue;
        
        double y = part.rapidity();
        double pt = part.pT();

        _pt_hist->fill(y, pt);
        if (y > 2.5 && y < 4.0)
          _forw_y_pt_hist->fill(pt);
        else if (y > -4.0 && y < -2.5)
          _backw_y_pt_hist->fill(pt);
      }
    }

    void finalize()
    {
      // Compute scale factor with inelastic cross-section from input file and sum of
      // weights (corresponds to number of events in input file)
      const double scale_fact = crossSection() / millibarn / sumOfWeights();

      // Apply scale factor (histogram scaled by transverse-momentum interval widths
      // when running plotting command [results in differential cross-section])
      _pt_hist->divByGroupWidth();
      _pt_hist->scaleW(scale_fact);

      // Compute forward-to-backward ratio
      divide(_forw_y_pt_hist, _backw_y_pt_hist, _forw_to_backw_rat);
    }

    Histo1DGroupPtr _pt_hist;
    Estimate1DPtr _forw_to_backw_rat;
    Histo1DPtr _forw_y_pt_hist;
    Histo1DPtr _backw_y_pt_hist;
  };

  RIVET_DECLARE_PLUGIN(LHCB_2022_I2694685);
}
