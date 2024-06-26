// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> → D* D* pi0
  class BESIII_2015_I1382146 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1382146);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      // histograms
      // histograms
      if (isCompatibleWithSqrtS(4.23)) {
        _ecms="4.23";
        book(_h[0],2,1,2);
      }
      else if (isCompatibleWithSqrtS(4.26)) {
        _ecms="4.26";
        book(_h[0],2,1,3);
      }
      else {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
      book(_h[1],2,1,1);
      book(_sigma,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DD,other;
      for (const Particle & p : fs) {
        Particle parent=p;
        while (!parent.parents().empty()) {
          parent=parent.parents()[0];
          if(parent.abspid()==413 || parent.abspid()==423) break;
        }
        if (parent.abspid()==413 || parent.abspid()==423) {
          bool found=false;
          for (const auto& D : DD) {
            // D already in list
            if (fuzzyEquals(D.mom(), parent.mom())) {
              found=true;
              break;
            }
          }
          if (!found) DD.push_back(parent);
        }
        else {
          other.push_back(p);
        }
      }
      // D Dbar + neutral pion
      if (DD.size()!=2 || other.size()!=1) vetoEvent;
      if (DD[0].pid()!=-DD[1].pid()) vetoEvent;
      if (other[0].pid()!=111) vetoEvent;
      const double mass = (DD[0].momentum()+DD[1].momentum()).mass();
      _h[0]->fill(mass/GeV);
      _h[1]->fill(mass/GeV);
      // parent Z0
      if (DD[0].parents()[0].pid()==_pid && DD[1].parents()[0].pid()==_pid &&
          fuzzyEquals(DD[0].parents()[0].momentum(),DD[1].parents()[0].momentum()) ) _sigma->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // distributions
      normalize(_h, 1.0, false);
      // cross section
      scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    Histo1DPtr _h[2];
    int _pid;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1382146);

}
