// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  /// @brief Measurement of b-hadron fractions in 13 TeV pp collisions
  class LHCB_2019_I1720859 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1720859);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Unstable particles
      const UnstableParticles up(Cuts::absetaIn(2, 5) && Cuts::ptIn(4.*GeV, 25.*GeV) && (Cuts::abspid == 511 || Cuts::abspid == 521 || Cuts::abspid == 531 || Cuts::abspid == 5122));
      declare(up, "up");

      // Book temp histograms
      book(_h["Bx"], "TMP/Bx", refData(1, 1, 1));
      book(_h["Bs"], "TMP/Bs", refData(1, 1, 1));
      book(_h["Lb"], "TMP/Lb", refData(1, 1, 1));
      // Bs and Lb ratios to B mesons
      book(_s["BsRatio"], 1, 1, 1);
      book(_s["LbRatio"], 1, 1, 2);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const UnstableParticles& up = apply<UnstableParticles>(event, "up");

      // Fill histograms with B hadron pTs
      for (const Particle& p : up.particles()) {
        if (p.abspid() == 521 || p.abspid() == 511)
          _h["Bx"]->fill(p.pT() / GeV);
        else if (p.abspid() == 531)
          _h["Bs"]->fill(p.pT() / GeV);
        else if (p.abspid() == 5122)
          _h["Lb"]->fill(p.pT() / GeV);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      divide(_h["Bs"], _h["Bx"], _s["BsRatio"]);
      divide(_h["Lb"], _h["Bx"], _s["LbRatio"]);
    }

    /// @}

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Estimate1DPtr> _s;
    /// @}
  };

  RIVET_DECLARE_PLUGIN(LHCB_2019_I1720859);
}
