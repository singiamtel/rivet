// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Kaon spectra at 3.63, 4.03 and 4.5 GeV
  class PLUTO_1977_I118873 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1977_I118873);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::pid==PID::K0S), "UFS");
      book(_h_sigma,1,1,1);
      for (const string& en :_h_sigma.binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms = en;
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms = en;
            break;
          }
        }
      }
      if (isCompatibleWithSqrtS(3.63*GeV)) {
        book(_h_spectrum, 2, 1, 1);
      }
      else if (isCompatibleWithSqrtS(4.03*GeV)) {
        book(_h_spectrum, 3, 1, 1);
      }
      else if (isCompatibleWithSqrtS(4.5*GeV)) {
        book(_h_spectrum, 4, 1, 1);
      }
      else if(_ecms.empty())
        MSG_ERROR("Beam energy not supported!");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = 0.5*(beams.first.p3().mod() + beams.second.p3().mod());
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // unstable particles
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        _h_sigma->fill(_ecms);
        const double xp = p.E()/meanBeamMom;
        const double beta = p.p3().mod()/p.E();
        if(_h_spectrum) _h_spectrum->fill(xp,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_spectrum)
        scale(_h_spectrum, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale(_h_sigma, crossSection()/nanobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spectrum;
    BinnedHistoPtr<string> _h_sigma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1977_I118873);


}
