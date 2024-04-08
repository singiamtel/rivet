// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief D*+ production at 34.4 GeV
  class JADE_1984_I202785 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1984_I202785);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_x, 1, 1, 1);
      book(_h_theta, 3, 1, 1);
      _axis = YODA::Axis<double>(6, -0.9, 0.9);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty())  _edges = _h_theta->xEdges();
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      Vector3 axis = beams.first.pid() == PID::EMINUS ? beams.first.p3().unit() : beams.second.p3().unit();
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==413)) {
        double modp = p.p3().mod();
        double xE =p.E()/meanBeamMom;
        _h_x->fill(xE);
        if(xE>0.4 && p.pid()>0) {
          _h_theta->fill(map2string(p.p3().dot(axis)/modp));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_theta,33.6,false); // normalize to data
      scale(_h_x, crossSection()/microbarn/sumOfWeights()*sqr(sqrtS()));
      for(auto & b: _h_theta->bins()) {
        const size_t idx = b.index();
        b.scaleW(1./_axis.width(idx));
      }

    }

    /// @}

    string map2string(const double value) const {
      const size_t idx = _axis.index(value);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }


    /// @name Histograms
    /// @{
    Histo1DPtr _h_x;
    BinnedHistoPtr<string> _h_theta;
    YODA::Axis<double> _axis;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(JADE_1984_I202785);


}
