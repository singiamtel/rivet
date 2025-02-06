// -*- C++ -*-
#include "Rivet/Analysis.hh"

namespace Rivet {

/// @brief An example demonstrating the usage of different jet collections and reclustering.
class EXAMPLE_JETS : public Analysis {
public:

  /// Constructor
  RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_JETS);


  /// @name Analysis methods
  /// @{

  /// Book histograms and initialise projections before the run
  void init() {
    FinalState fs(Cuts::abseta < 5);
    FastJets fj(fs, JetAlg::ANTIKT, 0.4);
    _smallR_jDef = fj.jetDef();
    declare(fj, "smallRJets");

    std::vector<std::string> jetTypes = {"smallRJets", "trimJets", "reclusterJets", "trackJets", "sisconeJets"};
    for ( const std::string &s : jetTypes ){
      book(_h["N"+s],  "N"+s,  10, 0, 10);
      book(_h["NB"+s], "NB"+s, 10, 0, 10);
      book(_h[s+"Pt"], s+"Pt", 20, 0, 1000);
    }
  }

  /// Perform the per-event analysis
  void analyze(const Event& event) {
    std::unordered_map<string, Jets> jets;

    // Small-R jets example
    jets["smallRJets"] = apply<JetFinder>(event, "smallRJets").jetsByPt(Cuts::pT > 10*GeV);

    // Jet trimming example
    jets["trimJets"] = trimJetsFrac(jets["smallRJets"], 0.1, cmpMomByMass); // default is sorting by pT

    // Reclustering jets example from enum, optionally with filter
    const fastjet::JetDefinition _largeR_jDef = FastJets::mkJetDef(JetAlg::ANTIKT, 0.8);
    const fastjet::Filter filter(_smallR_jDef /*or _largeR_jDef*/, fastjet::SelectorPtFractionMin(0.1));
    jets["reclusterJets"] = reclusterJets<JetAlg::VARIABLER>(jets["smallRJets"], /*filter,*/ 550*GeV, 0.4, 1.5, fastjet::contrib::VariableRPlugin::AKTLIKE);
    std::vector<uint> nPJConstsVec;
    for( const Jet &jet : jets["reclusterJets"] ){
      // Jets in reclusterJets() will have the original jet constituents (=hadrons) as constituents
      // Access the number of jets that have been reclustered into a new variable-R jet like this:
      uint nPJConsts = jet.pseudojet().constituents().size();
      nPJConstsVec.push_back(nPJConsts);
    }

    // Reclustering jets example from object (enum as arg, radius as arg), optionally with filter
    const FJPluginPtr trackjetPlugin = FastJets::mkPlugin(JetAlg::TRACKJET, 0.4);
    jets["trackJets"] = reclusterJets(jets["smallRJets"], trackjetPlugin/*, filter*/); // already sorted by pT as part of reclustering

    // Reclustering jets example from object (enum as template, arbitrary args), optionally with filter
    const FJPluginPtr sisconePlugin = FastJets::mkPlugin<JetAlg::SISCONE>(0.4, 0.75);
    jets["sisconeJets"] = reclusterJets(jets["smallRJets"], sisconePlugin/*, filter*/);

    // Retrieve b-jet multiplicities
    map<string, size_t> bJetsMult;
    for( auto const& [key, val] : jets ){
      const Jets bJets = select(val, hasBTag(Cuts::pT>5*GeV, 0.3));
      bJetsMult[key] = bJets.size();
    }

    // Fill histograms
    for( auto const& [key, val] : jets ){
      _h["N"+key]->fill(val.size());
      _h["NB"+key]->fill(bJetsMult[key]);
      for( const Jet &jet : val ){
        _h[key+"Pt"]->fill(jet.pT());
      }
    }
  }


  /// Normalise histograms etc., after the run
  void finalize() {
    scale(_h, crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
  }

  /// @}


  private:

  /// @name Aux objects
  /// @{
  fastjet::JetDefinition _smallR_jDef;
  /// @}

  /// @name Histograms
  /// @{
  map<string, Histo1DPtr> _h;
  /// @}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(EXAMPLE_JETS);
}
