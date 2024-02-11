// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// Measurement of the W and Z bosons pT produced in pp collisions at sqrt(s)=8 TeV
  class CMS_2016_I1471281 : public Analysis {
  public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2016_I1471281)


      /// @name Analysis methods
      /// @{

      /// Book histograms and initialise projections before the run
      void init() {

        // Get options from the new option system
        // default to both.
        if ( getOption("VMODE") == "BOTH" ) _mode = 0;
        if ( getOption("VMODE") == "W" )    _mode = 1;
        if ( getOption("VMODE") == "Z" )    _mode = 2;

        // Set up projections
        Cut cut_mu = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;

        // Dressed W's
        PromptFinalState mufs(cut_mu && Cuts::abspid == PID::MUON);
        declare(mufs, "Muons");
        declare(MissingMomentum(), "MET");

        // Dressed Z's
	/// @note There's no protection in this analysis to ensure the Z and W don't use the same muon
        DileptonFinder zmumu_Finder(91.2*GeV, 0.0, cut_mu && Cuts::abspid == PID::MUON, Cuts::massIn(60*GeV, 120*GeV)); //< no dressing
        declare(zmumu_Finder, "Zmumu_Finder");

        // Histograms
        if (_mode == 0 || _mode == 1) {
          book(_hist_WtoMuNuPt, 8, 1, 1);
        }
        if (_mode == 0 || _mode == 2) {
          book(_hist_ZtoMuMuPt, 9, 1, 1);
        }
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        if (_mode == 0 || _mode == 1) {
          // Get the W bosons - muon decay channel
          const P4& pmiss = apply<MissingMomentum>(event, "MET").missingMom();
          const Particles& mus = apply<PromptFinalState>(event, "Muons").particles();
          const Particles mus_mtfilt = select(mus, [&](const Particle& m){ return mT(m, pmiss) > 0*GeV; });
          const int imfound = closestMatchIndex(mus_mtfilt, pmiss, Kin::mass, 80.4*GeV);
          if (imfound >= 0) {
            const FourMomentum pWmunu = mus_mtfilt[imfound] + pmiss;
            _hist_WtoMuNuPt->fill(pWmunu.pT()/GeV);
          }
        }

        if (_mode == 0 || _mode == 2) {
          // Get the Z bosons - muon decay channel
	  /// @todo Note overlap with the W muon in diboson events
          const DileptonFinder& zmumu_Finder = apply<DileptonFinder>(event, "Zmumu_Finder");
          if (!zmumu_Finder.bosons().empty()) {
            const FourMomentum pZmumu = zmumu_Finder.bosons()[0].momentum();
            _hist_ZtoMuMuPt->fill(pZmumu.pT()/GeV);
          }
        }

      }


      /// Normalise histograms etc., after the run
      void finalize() {
        MSG_DEBUG("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection()/picobarn << " pb");
        MSG_DEBUG("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
        MSG_DEBUG("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
        if (_mode == 0 || _mode == 1) normalize(_hist_WtoMuNuPt);
        if (_mode == 0 || _mode == 2) normalize(_hist_ZtoMuMuPt);
      }

      /// @}


  private:

      // Data members
      size_t _mode = 0;
    

      /// @name Histograms
      /// @{
      Histo1DPtr _hist_WtoMuNuPt;
      Histo1DPtr _hist_ZtoMuMuPt;
      /// @}
    
  };


  RIVET_DECLARE_PLUGIN(CMS_2016_I1471281);

}
