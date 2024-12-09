// -*- C++ -*-
#include "Rivet/SimpleAnalysis.hh"

namespace Rivet {


  /// @brief Demonstrate the SimpleAnalysis interface
  class EXAMPLE_SIMPLEANALYSIS : public SimpleAnalysis {
  public:

    /// Constructor
    RIVET_DEFAULT_SIMPLEANALYSIS_CTOR(EXAMPLE_SIMPLEANALYSIS);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {

      // Set the detector-smearing family
      setDetSmearing(ATLAS_RUN2, MV2C20);

      // Book histograms
      book(_h["ne"], "num_elec", 7, -0.5, 6.5);
      book(_h["nm"], "num_muon", 7, -0.5, 6.5);
      book(_h["nt"], "num_tau", 7, -0.5, 6.5);
      book(_h["ny"], "num_photon", 7, -0.5, 6.5);
      book(_h["njb"], "num_bjet", 7, -0.5, 6.5);
      book(_h["ptjb1"], "pt_bjet1", 20, 30, 200);
      book(_h["ptl1"], "pt_lep1", 20, 20, 200);
      // Book counter
      book(_c["npass"], "count_pass");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Apply a missing-momentum cut
      if (met().mod() < 100*GeV) vetoEvent;

      // Run the overlap removal (see the SimpleAnalysis base class for some info)
      doSimpleOverlapRemoval();

      // Debug writing
      MSG_DEBUG("Numbers of e, mu, tau, photon,and b-jet: "
                << electrons().size() << ", "
                << muons().size() << ", "
                << taus().size() << ", "
                << photons().size() << ", "
                << bjets().size());

      // Veto event if there are no b-jets (example eta cut to fully contain a 0.4 jet in the tracker)
      const Jets jbs = bjets(Cuts::pT > 30*GeV && Cuts::abseta < 2.1);
      if (jbs.empty()) vetoEvent;

      // Get electrons or muons and require one to three of them
      const Particles leps =
        sortByPt(electrons(Cuts::pT > 20*GeV && Cuts::abseta < 2.5) +
                 muons(Cuts::pT > 20*GeV && Cuts::abseta < 2.5));

      // Fill histograms with object counts, and leading b-jet and lepton pTs
      _h["ne"]->fill(electrons().size());
      _h["nm"]->fill(muons().size());
      _h["nt"]->fill(taus().size());
      _h["ny"]->fill(photons().size());
      _h["njb"]->fill(jbs.size());
      if (!jbs.empty()) _h["ptjb1"]->fill(jbs[0].pT()/GeV);
      if (!leps.empty()) _h["ptl1"]->fill(leps[0].pT()/GeV);

      // Fill counter for (made-up) signal-region passes
      if (leps.size() != 1 && bjets(Cuts::pT > 100*GeV).size() == 2) {
        _c["npass"]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h); // normalize to unity
      scaleToIntLumi(_c, 300/femtobarn); // normalize to 300/fb
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_SIMPLEANALYSIS);

}
