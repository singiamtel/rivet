// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief H1 charged particle photoproduction
  class H1_1998_I477556 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1998_I477556);

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      // Acceptance
      const FinalState fs(Cuts::abseta < 1 && Cuts::pT > 2*GeV);
      const ChargedFinalState cfs(fs);
      declare(cfs, "CFS");

      // Table 1, pT spectra
      book(_h_ch_pT, 1, 1, 1);
      // Table 2, eta distributions with different pT cuts.
      book(_h_ch_eta[0], 2, 1, 1);
      book(_h_ch_eta[1], 2, 1, 2);
    }

    // Do the analysis
    void analyze(const Event& event) {

      // Loop over charged particles and fill histograms.
      const Particles& particles = apply<FinalState>(event, "CFS").particles();
      for (const Particle& p : particles) {

        // Cross section for dsigma/dpt^2deta but as a function of pT.
        // Use the conversion factor as a weight.
        _h_ch_pT->fill(p.pt()/GeV, 1.0 / (2*p.pt()/GeV));
        _h_ch_eta[0]->fill(p.eta());
        if (p.pt() > 3. * GeV) _h_ch_eta[1]->fill(p.eta());
      }
    }

    // Finalize
    void finalize() {
      const double sf = crossSection()/nanobarn/sumOfWeights();
      // For d^2sigma/dpT^2deta scale with eta bin width.
      scale(_h_ch_pT, 0.5 * sf);
      scale(_h_ch_eta, sf);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ch_pT, _h_ch_eta[2];
    //@}

  };

  RIVET_DECLARE_PLUGIN(H1_1998_I477556);

}
