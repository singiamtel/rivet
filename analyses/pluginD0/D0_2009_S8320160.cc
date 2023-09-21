// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 dijet angular distributions
  class D0_2009_S8320160 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2009_S8320160);


    /// @name Analysis methods
    /// @{

    // Book histograms
    void init() {
      FinalState fs;
      FastJets conefinder(fs, FastJets::D0ILCONE, 0.7);
      declare(conefinder, "ConeFinder");

      book(_h_chi_dijet, {250., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1960.});
      for (auto& b : _h_chi_dijet->bins()) {
        book(b, b.index(), 1, 1);
      }
    }


    /// Do the analysis
    void analyze(const Event & e) {

      const Jets& jets = apply<JetAlg>(e, "ConeFinder").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();

      if (fabs(y0+y1)>2) vetoEvent;

      double mjj = FourMomentum(j0+j1).mass();
      double chi = exp(fabs(y0-y1));
      if(chi<16.)  _h_chi_dijet->fill(mjj, chi);
    }


    /// Finalize
    void finalize() {
      normalize(_h_chi_dijet);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_chi_dijet;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2009_S8320160, D0_2009_I824127);

}
