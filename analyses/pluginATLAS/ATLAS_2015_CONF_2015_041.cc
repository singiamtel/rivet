// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Z + jets in pp at 13 TeV
  ///
  /// @note This base class contains a "mode" variable for combined, e, and mu channel derived classes
  class ATLAS_2015_CONF_2015_041 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2015_CONF_2015_041);

    /// @}


    /// Book histograms and initialise projections before the run
    void init() {


      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      else if ( getOption("LMODE") == "MU" ) _mode = 2;

      Cut cuts = Cuts::pT > 25*GeV && Cuts::abseta < 2.5;
      DileptonFinder eefinder(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::ELECTRON,
                             Cuts::massIn(66*GeV, 116*GeV));
      DileptonFinder mmfinder(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::MUON,
                             Cuts::massIn(66*GeV, 116*GeV));
      declare(eefinder, "eefinder");
      declare(mmfinder, "mmfinder");

      // Define veto FS in order to prevent Z-decay products entering the jet algorithm
      VetoedFinalState had_fs;
      had_fs.addVetoOnThisFinalState(eefinder);
      had_fs.addVetoOnThisFinalState(mmfinder);
      FastJets jets(had_fs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jets, "jets");

      // combination
      book(_h["Njets_comb"], 1, 1, 1);
      book(_e["ratio_comb"], 2, 1, 1);
      // individual channels
      if (_mode == 0 || _mode == 1) {
        book(_h["Njets_el"], 1, 2, 1);
        book(_e["ratio_el"], 2, 2, 1);
      }
      if (_mode == 0 || _mode == 2) {
        book(_h["Njets_mu"], 1, 2, 2);
        book(_e["ratio_mu"], 2, 2, 2);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& elecs = apply<DileptonFinder>(event, "eefinder").constituents();
      const Particles& muons = apply<DileptonFinder>(event, "mmfinder").constituents();

      if (elecs.size() + muons.size() != 2)  vetoEvent;

      Jets jets = apply<JetFinder>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 2.5);
      idiscardIfAnyDeltaRLess(jets, elecs, 0.4);
      idiscardIfAnyDeltaRLess(jets, muons, 0.4);

      fillHistos(elecs, jets);
      fillHistos(muons, jets);
    }

    void fillHistos(const Particles& leptons, const Jets& jets) {
      if (leptons.size() != 2)  return;

      const size_t njets = jets.size();
      for (size_t i = 0; i <= njets; ++i) {
        if (_mode == 0 || _mode == 1) {
          _h["Njets_comb"]->fill(i + 0.5);
          _h["Njets_el"]->fill(i + 0.5);
        }
        if (_mode == 0 || _mode == 2) {
          _h["Njets_comb"]->fill(i + 0.5);
          _h["Njets_mu"]->fill(i + 0.5);
        }
      }

    }

    void finalize() {
      scale(_h, crossSectionPerEvent());
      if (_mode == 0)  scale(_h["Njets_comb"], 0.5); // average of el + mu

      vector<string> channels = { "comb" };
      if (_mode == 0 || _mode == 1)  channels.push_back("el");
      if (_mode == 0 || _mode == 2)  channels.push_back("mu");
      YODA::Estimate num, den;
      for (size_t i = 1; i < 5; ++i) {
        for (const string& ch : channels) {
          num.set(_h["Njets_"+ch]->bin(i+1).sumW(), _h["Njets_"+ch]->bin(i+1).errW());
          den.set(_h["Njets_"+ch]->bin(i).sumW(),   _h["Njets_"+ch]->bin(i).errW());
          _e["ratio_"+ch]->bin(i) = YODA::divide(num, den);
        }
      }

    }

    /// @}


  private:

    size_t _mode;

    map<string, Estimate1DPtr> _e;
    map<string, Histo1DPtr> _h;
  };


  RIVET_DECLARE_PLUGIN(ATLAS_2015_CONF_2015_041);
}
