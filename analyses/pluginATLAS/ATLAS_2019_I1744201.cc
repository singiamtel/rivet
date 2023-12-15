// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet{


  /// Z+jet at 8 TeV
  class ATLAS_2019_I1744201 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1744201);


    void init() {

      DileptonFinder zfinder_el(91.2*GeV, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 20*GeV &&
                                Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
      declare(zfinder_el, "DileptonFinder_el");

      declare(FastJets(zfinder_el.remainingFinalState(), JetAlg::ANTIKT, 0.4,
                       JetMuons::NONE, JetInvisibles::NONE), "AKT04");

      h_jet_y_pt.resize(6);
      for (size_t iPtBin=0; iPtBin < h_jet_y_pt.size(); ++iPtBin) {
        book(h_jet_y_pt[iPtBin], iPtBin+2, 1, 1);
      }

    }


    void analyze(const Event& event) {

      // electrons selection
      const DileptonFinder& zfinder = apply<DileptonFinder>(event, "DileptonFinder_el");
      if ( zfinder.bosons().size() != 1)  vetoEvent;

      const Particles& leptons = zfinder.constituents();
      if ( leptons.size() != 2)  vetoEvent;

      if (deltaR(leptons[0], leptons[1]) < 0.2)  vetoEvent;


      // jets selection
      Jets jets = apply<FastJets>(event, "AKT04").jetsByPt( Cuts::pT > 25*GeV && Cuts::absrap < 3.4 );
      idiscardIfAnyDeltaRLess(jets, leptons, 0.4);
      if (jets.empty())  vetoEvent;  // require at least one jet in event

      for (const Jet& jet : jets) {
        const double jet_pt = jet.pT() / GeV;
        for(size_t iPtBin = 0; iPtBin < (ptBins.size() - 1); ++iPtBin) {
          if (jet_pt >= ptBins[iPtBin] && jet_pt < ptBins[iPtBin+1]) {
            h_jet_y_pt[iPtBin]->fill(jet.absrap());
          }
        }
	    }
   }

  void finalize() {

    const double norm = crossSection()/femtobarn/sumOfWeights();
    for(int iPtBin=0; iPtBin < 6; ++iPtBin){
    	scale(h_jet_y_pt[iPtBin], norm / ( ptBins[iPtBin+1] - ptBins[iPtBin] ));
    }

  }

   protected:

     vector<double> ptBins = { 25., 50., 100., 200., 300., 400., 1050. };

   private:
     vector<Histo1DPtr> h_jet_y_pt;

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1744201);
}
