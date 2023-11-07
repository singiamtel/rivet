// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Inclusive jets at 8 TeV
  class ATLAS_2017_I1604271 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2017_I1604271);

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs,"FinalState");
      FastJets fj04(fs, JetAlg::ANTIKT, 0.4);
      FastJets fj06(fs, JetAlg::ANTIKT, 0.6);
      fj04.useInvisibles();
      declare(fj04, "AntiKT04");
      fj06.useInvisibles();
      declare(fj06, "AntiKT06");

      // |y| and ystar bins
      const vector<double> ybins{ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

      // Book histograms
      // pT histograms
      book(_pThistograms6, ybins);
      book(_pThistograms4, ybins);
      for (size_t i = 0; i < _pThistograms4->numBins(); ++i) {
        book(_pThistograms6->bin(i+1), i+1, 1, 1);
        book(_pThistograms4->bin(i+1), i+7, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& kt4Jets = apply<FastJets>(event, "AntiKT04").jetsByPt(Cuts::pT>70*GeV && Cuts::absrap < 3.0);
      const Jets& kt6Jets = apply<FastJets>(event, "AntiKT06").jetsByPt(Cuts::pT>70*GeV && Cuts::absrap < 3.0);

      int nJets4 = kt4Jets.size();
      int nJets6 = kt6Jets.size();

      // Inclusive jet selection
      for(int ijet=0;ijet<nJets4;++ijet){ // loop over jets
        FourMomentum jet = kt4Jets[ijet].momentum();
        // pT selection
        const double absy = jet.absrap();
        _pThistograms4->fill(absy,jet.pt()/GeV);
      }

      for(int ijet=0;ijet<nJets6;++ijet){ // loop over jets
        FourMomentum jet = kt6Jets[ijet].momentum();
        // pT selection
        const double absy = jet.absrap();
        _pThistograms6->fill(absy,jet.pt()/GeV);
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs_norm_factor = 0.5*crossSection() / picobarn / sumOfWeights();

      scale(_pThistograms4, xs_norm_factor);
      scale(_pThistograms6, xs_norm_factor);
      divByGroupWidth({_pThistograms4, _pThistograms6});

    }

  private:

    // The inclusive pT spectrum for akt4 jets
    Histo1DGroupPtr _pThistograms4;
    // The inclusive pT spectrum for akt6 jets
    Histo1DGroupPtr _pThistograms6;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2017_I1604271);

}
