// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Rdphi observable for the determination of the strong coupling at 13 TeV
  class CMS_2024_I2780732 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2024_I2780732);

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs;
      FastJets akt(fs, JetAlg::ANTIKT, 0.7);
      declare(akt, "antikT");

      // Book histograms
      // Denominator and numerator of Rdphi observable
      book(_h_Inclusive_Denominator, "_Inclusive_Denominator", refData(2, 1, 1));
      book(_h_Dphi_Numerator, "_Inclusive_Numerator", refData(2, 1, 1));
      book(_h_Rdphi, 2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Select jets with pt > 50 GeV and absolute rapidity < 2.5
      const Jets& jets = apply<FastJets>(event, "antikT").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 50.*GeV);

      // Loop over all jets and fill the denominator (inclusive jets)
      for (unsigned int i=0; i<jets.size(); ++i) {
        _h_Inclusive_Denominator->fill(jets[i].pT()/GeV);

        // Loop looking for neighboring jets for each jet
        for (unsigned int j=0; j<jets.size(); ++j) {
          if (i == j) continue;
          const double dphi = deltaPhi(jets[i], jets[j]);

          // Select neighboring jets and fill the numerator
          if ( (dphi>2.*pi/3.) && (dphi<7.*pi/8.) && (jets[j].pT() > 100.*GeV) ) {
            _h_Dphi_Numerator->fill(jets[i].pT()/GeV);
          }
        }
      }

    }


    /// Calculation of Rdphi observable
    void finalize() {

      divide(_h_Dphi_Numerator, _h_Inclusive_Denominator, _h_Rdphi);

    }


    /// @name Histograms
    Histo1DPtr _h_Dphi_Numerator, _h_Inclusive_Denominator;
    Estimate1DPtr _h_Rdphi;


  };


  RIVET_DECLARE_PLUGIN(CMS_2024_I2780732);

}
