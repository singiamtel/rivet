// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  kaon production at low energies
  class BESIII_2021_I1868813 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1868813);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_c_kaons   , "TMP/nK", refData<YODA::BinnedEstimate<string>>(1, 1, 1));

      for (const string& en : _c_kaons.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      unsigned int count = ufs.particles(Cuts::pid==310).size();
      _c_kaons->fill(_ecms,count);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> mult;
      book(mult,1, 1, 1);
      barchart(_c_kaons,mult);
      scale(mult,crossSection()/nanobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedProfilePtr<string> _c_kaons;
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1868813);

}
