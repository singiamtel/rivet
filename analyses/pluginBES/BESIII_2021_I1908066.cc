// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi production at low energies
  class BESIII_2021_I1908066 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1908066);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_c_psi   , 1, 1, 1);

      for (const string& en : _c_psi.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      unsigned int count = ufs.particles(Cuts::pid==443).size();
      _c_psi->fill(_ecms,count);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c_psi, crossSection()/nanobarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _c_psi;
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1908066);

}
