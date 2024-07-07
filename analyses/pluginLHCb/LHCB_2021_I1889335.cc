// -*-C++ - *-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {


  /// @brief Charged particle production at 13 TeV
  class LHCB_2021_I1889335 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2021_I1889335);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Register projection for primary particles
      declare(ALICE::PrimaryParticles(Cuts::etaIn(2.0, 4.8) && Cuts::abscharge > 0), "APRIM");

      vector<double> edges = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 4.8};
      book(_h_ppInel_neg, edges);
      book(_h_ppInel_pos, edges);
      for (size_t i=1; i < edges.size(); ++i) {
        book(_h_ppInel_neg->bin(i), 1, 1, i);
        book(_h_ppInel_pos->bin(i), 2, 1, i);
      }
    }


    void analyze(const Event &event) {

      const Particles cfs = apply<ALICE::PrimaryParticles>(event, "APRIM").particles();

      for (const Particle& myp : cfs) {
        if (myp.charge() < 0) {
          _h_ppInel_neg->fill(myp.eta(), myp.pT()/GeV);
        }
        else {
          _h_ppInel_pos->fill(myp.eta(), myp.pT()/GeV);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double scale_factor = crossSection() / millibarn / sumOfWeights();
      scale({_h_ppInel_neg, _h_ppInel_pos}, scale_factor);
      divByGroupWidth({_h_ppInel_neg, _h_ppInel_pos});
    }

    /// @}


  private:

    /// @name Histogram
    Histo1DGroupPtr _h_ppInel_neg, _h_ppInel_pos;

  };


  RIVET_DECLARE_PLUGIN(LHCB_2021_I1889335);

}
