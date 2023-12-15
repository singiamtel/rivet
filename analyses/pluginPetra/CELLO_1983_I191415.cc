// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi0 and gamma spectra at 14, 22 and 34
  class CELLO_1983_I191415 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1983_I191415);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      unsigned int iloc(0);
      if (isCompatibleWithSqrtS(14*GeV)) {
        iloc=1;
        _axes["gamma"] = YODA::Axis<double>({0.0428, 0.057, 0.07125, 0.0855, 0.09975, 0.114, 0.12825,
                                             0.14255, 0.17105, 0.19955, 0.22805, 0.25655, 0.2855,
                                             0.35625, 0.42755, 0.4988, 0.57, 0.8552});
        _axes["pi0"] = YODA::Axis<double>({0.1325, 0.1635, 0.2505, 0.3375, 0.434, 0.54, 0.791, 1.04});
      }
      else if (isCompatibleWithSqrtS(22*GeV)) {
        iloc=2;
        _axes["gamma"] = YODA::Axis<double>({0.02725, 0.03635, 0.04545, 0.05455, 0.06365, 0.07275,
                                             0.08185, 0.09095, 0.1091, 0.1273, 0.1455, 0.1637, 0.1819,
                                             0.22735, 0.27285, 0.31835, 0.36385, 0.54575, 0.72765});
        _axes["pi0"] = YODA::Axis<double>({0.0825, 0.1055, 0.1605, 0.2155, 0.2705, 0.3555, 0.5145, 0.733, 0.961});
      }
      else if (isCompatibleWithSqrtS(34*GeV)) {
        iloc=3;
        _axes["gamma"] = YODA::Axis<double>({0.01765, 0.02355, 0.02945, 0.0353, 0.04115, 0.04705, 0.05295,
                                             0.05885, 0.0706, 0.08235, 0.0941, 0.106, 0.1177, 0.1471, 0.1765,
                                             0.2059, 0.2353, 0.35295, 0.47065});
        _axes["pi0"] = YODA::Axis<double>({0.054, 0.068, 0.104, 0.14, 0.176, 0.23, 0.338, 0.45, 0.71, 0.97});
      }
      else {
        MSG_ERROR("Beam energy not supported!");
      }


      // Book histograms
      book(_h["gamma"], iloc,   1, 1);
      book(_h["pi0"],   iloc+3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        for (const auto& item : _h) {
          _edges[item.first] = item.second->xEdges();
        }
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
      const double meanBeamMom = 0.5*(beams.first.p3().mod() + beams.second.p3().mod());
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state to get particle spectra
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==111)) {
        double xE = p.E()/meanBeamMom;
        discfill("pi0", xE);
      }
      for (const Particle& p : apply<FinalState>(event, "FS").particles(Cuts::pid==22)) {
        double xE = p.E()/meanBeamMom;
        discfill("gamma", xE);
      }
    }

    void discfill(const string& name, const double value) {
      string edge = "OTHER";
      const size_t idx = _axes[name].index(value);
      if (idx && idx <= _edges[name].size())  edge = _edges[name][idx-1];
      _h[name]->fill(edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double fact = sqr(sqrtS())/GeV2*crossSection()/microbarn/sumOfWeights();
      scale(_h, fact);

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, BinnedHistoPtr<string>> _h;
    map<string, YODA::Axis<double>> _axes;
    map<string, vector<string>> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CELLO_1983_I191415);


}
