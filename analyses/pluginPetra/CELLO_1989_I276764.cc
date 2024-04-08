// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CELLO_1989_I276764 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1989_I276764);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==111 || Cuts::pid==221), "UFS");
      // Book histograms
      for(unsigned int ix=0;ix<4;++ix)
        book(_h[ix], 2+ix, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_edges[0].empty()) {
        for(unsigned int ix=0;ix<4;++ix)
          _edges[ix] = _h[ix]->xEdges();
      }
      // at least 5 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      if (numParticles < 5) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state to get particle spectra
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        double xE = p.E()/meanBeamMom;
        if(p.pid()==111) {
          _h[1]->fill( map2string(xE,1));
          _h[2]->fill( map2string(xE,2));
        }
        else
          _h[3]->fill( map2string(xE,3));
      }
      for (const Particle& p : apply<FinalState>(event, "FS").particles(Cuts::pid==22)) {
        _h[0]->fill(map2string(p.E()/meanBeamMom,0));
      }
    }

    string map2string(const double val, const size_t axis) const {
      const size_t idx = _axes[axis].index(val);
      if (idx && idx <= _edges[axis].size())  return _edges[axis][idx-1];
      return "OTHER";
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = sqr(sqrtS())/GeV2*crossSection()/microbarn/sumOfWeights();
      for(unsigned int ix=0;ix<4;++ix) {
        scale(_h[ix], fact);
        for(auto & b: _h[ix]->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axes[ix].width(idx));
        }
      }
    }
    
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h[4];
    vector<string> _edges[4];
    YODA::Axis<double> _axes[4] = { YODA::Axis<double>({0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,
                                                            0.175,0.195,0.215,0.235,0.255,0.275,0.295,0.315,0.3535,0.394,0.434,0.4735,
                                                            0.514,0.555,0.5965,0.641,0.691,0.74,0.785,0.839,0.901}),
                                    YODA::Axis<double>({0.039,0.059,0.079,0.109,0.139,0.193,0.277,0.439}),
                                    YODA::Axis<double>({0.2705,0.3415,0.412,0.482,0.552,0.6225,0.6935,0.7805,0.876,0.962}),
                                    YODA::Axis<double>({0.0225,0.0775,0.1325,0.2875,0.4625,0.9175}) };
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CELLO_1989_I276764);


}
