// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 dijet invariant mass measurement
  class D0_2010_I846483 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_I846483);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      FastJets conefinder(fs, JetAlg::D0ILCONE, 0.7);
      declare(conefinder, "ConeFinder");

      book(_h_m_dijet, {0., 0.4, 0.8, 1.2, 1.6, 2., 2.4});
      for (auto& b : _h_m_dijet->bins()) {
        book(b, b.index(), 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      const Jets& jets = apply<JetFinder>(e, "ConeFinder").jetsByPt(Cuts::pT > 40.0*GeV);
      if (jets.size() < 2) vetoEvent;

      const double ymax = std::max(jets[0].absrap(), jets[1].absrap());
      const double mjj = FourMomentum(jets[0].mom() + jets[1].mom()).mass();

      _h_m_dijet->fill(ymax, mjj/TeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_m_dijet, crossSection()/sumOfWeights());
      divByGroupWidth(_h_m_dijet);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_m_dijet;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_I846483, D0_2010_S8566488);

}
