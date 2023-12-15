// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 measurement of the ratio \f$ \sigma(Z/\gamma^* + n \text{ jets})/\sigma(Z/\gamma^*) \f$
  class D0_2008_S6879055 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2008_S6879055);


    /// @name Analysis methods
    /// @{

    // Book histograms
    void init() {
      DileptonFinder zfinder(91.2*GeV, 0.2, Cuts::abspid == PID::ELECTRON, Cuts::massIn(40*GeV, 200*GeV));
      declare(zfinder, "DileptonFinder");

      FastJets conefinder(zfinder.remainingFinalState(), JetAlg::D0ILCONE, 0.5);
      declare(conefinder, "ConeFinder");

      book(_crossSectionRatio, 1, 1, 1);
      book(_pTjet1, 2, 1, 1);
      book(_pTjet2, 3, 1, 1);
      book(_pTjet3, 4, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {

      if (_edges.empty())  _edges = _crossSectionRatio->xEdges();

      const DileptonFinder& zfinder = apply<DileptonFinder>(event, "DileptonFinder");
      if (zfinder.bosons().size() != 1)  vetoEvent;

      FourMomentum e0 = zfinder.constituents()[0].mom();
      FourMomentum e1 = zfinder.constituents()[1].mom();

      Jets jets = apply<JetFinder>(event, "ConeFinder").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.5);
      idiscard(jets, deltaRLess(e0, 0.4));
      idiscard(jets, deltaRLess(e1, 0.4));

      // For normalisation of crossSection data (includes events with no jets passing cuts)
      _crossSectionRatio->fill(string("INC"));

      // Fill jet pT and multiplicities
      if (jets.size() >= 1) {
        _crossSectionRatio->fill(_edges[0]);
        _pTjet1->fill(jets[0].pT());
      }
      if (jets.size() >= 2) {
        _crossSectionRatio->fill(_edges[1]);
        _pTjet2->fill(jets[1].pT());
      }
      if (jets.size() >= 3) {
        _crossSectionRatio->fill(_edges[2]);
        _pTjet3->fill(jets[2].pT());
      }
      if (jets.size() >= 4) {
        _crossSectionRatio->fill(_edges[3]);
      }
    }


    /// Finalize
    void finalize() {
      // Now divide by the inclusive result (in the otherflow)
      scale(_crossSectionRatio, 1.0/_crossSectionRatio->bin(0).sumW());

      // Normalise jet pTs to integrals of data
      // @note There is no other way to do this, because these quantities are not detector-corrected
      normalize(_pTjet1, refData(2, 1, 1).auc());
      normalize(_pTjet2, refData(3, 1, 1).auc());
      normalize(_pTjet3, refData(4, 1, 1).auc());
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _crossSectionRatio;
    Histo1DPtr _pTjet1, _pTjet2, _pTjet3;
    vector<string> _edges;
    /// @}

  };


  RIVET_DECLARE_ALIASED_PLUGIN(D0_2008_S6879055, D0_2008_I724239);

}
