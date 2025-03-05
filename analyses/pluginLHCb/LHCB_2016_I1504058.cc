// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet
{

  /// @brief b-quark production in 7 and 13 TeV pp collisions 
  class LHCB_2016_I1504058 : public Analysis
  {

  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1504058);

    void init()
    {
      // Check energy of input file and book corresponding histogram
      if (isCompatibleWithSqrtS(7000 * GeV))
        en_case = 7;
      else if (isCompatibleWithSqrtS(13000 * GeV))
        en_case = 13;
      if (en_case == 0)
        throw UserError("Centre-of-mass energy of input file"
                        " is neither 7000 nor 13000 GeV.");
      book(_b_quark_hist[0], 1, 1, 1);
      book(_b_quark_hist[1], 2, 1, 1);
      // Select pseudorapidity range
      declare(UnstableParticles(Cuts::absetaIn(2.0, 5.0)), "UPs");
    }

    void analyze(const Event &ev)
    {
      const int hix = en_case == 7 ? 0 : 1;
      // Apply UnstableParticles projection to get entire decay chain of every particle
      const UnstableParticles &unst_parts = apply<UnstableParticles>(ev, "UPs");
      for (const Particle &part : unst_parts.particles())
      {
        // Select beauty hadrons
        if (!(PID::isBottomHadron(part.pid())))
          continue;

        // Skip bottomonia (excluded from definition) and Bc+ meson (ignored in
        // measurement)
        if ((PID::isQuarkonium(part.pid())) || (part.abspid() == PID::BCPLUS))
          continue;

        // Veto beauty hadrons w/ beauty hadron child (pick only children)
        if (part.isLastWith(hasBottom)) continue;

        // Double sample size by taking absolute value and then scale by 0.5 due to LHCb
        // being one-sided
        double eta = part.abseta();

        // Scale again by 0.5 to get mean of beauty quarks and antiquarks
        _b_quark_hist[hix]->fill(eta, 0.25);
      }
    }

    void finalize()
    {
      // Compute scale factor with inelastic cross-section from input file and sum of
      // weights (corresponds to number of events in input file) and do not modify it
      // for using symmetric LHCb phase space as taken into account by histogram count
      const double scale_fact = crossSection() / microbarn / sumOfWeights();

      // compute ratio of b quark productions at 13 and 7 TeV (during merge)
      if ((_b_quark_hist[0]->effNumEntries() > 0) &&  (_b_quark_hist[1]->effNumEntries() > 0)) 
      {
        Estimate1DPtr tmp;
        book(tmp, 3, 1, 1);
        divide(_b_quark_hist[1], _b_quark_hist[0], tmp);
      }

      // Apply scale factor (histogram scaled by interval widths when running plotting
      // command [results in differential cross-section])
      for (int i = 0; i < 2; i++) scale(_b_quark_hist[i], scale_fact);

    }

    Histo1DPtr _b_quark_hist[2];
    int en_case = -1;
  };

  RIVET_DECLARE_PLUGIN(LHCB_2016_I1504058);
}
