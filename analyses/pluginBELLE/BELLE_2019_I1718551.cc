// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, kaon and proton spectra at 10.58 GeV
  class BELLE_2019_I1718551 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1718551);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs, "FS");
      declare(Thrust(fs),"Thrust");
      declare(Beam(), "Beams");

      for (size_t ix=0; ix<6; ++ix) {

        if (ix == 1 || ix == 2) {
          book(_proton[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                             0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8});
        }
        else if (ix) {
          book(_proton[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                             0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85});
        }
        else {
          book(_proton[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65});
        }

        if (ix) {
          book(_pion[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                           0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0});
          book(_kaon[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                           0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0});
        }
        else {
          book(_pion[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                           0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85});
          book(_kaon[ix], {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                           0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85});
        }
        for (size_t iy=0; iy < _pion[ix]->numBins(); ++iy) {
          book(_pion[ix]->bin(iy+1), 1, ix+1, iy+1);
          book(_kaon[ix]->bin(iy+1), 2, ix+1, iy+1);
        }
        for (size_t iy=0; iy < _proton[ix]->numBins(); ++iy) {
          book(_proton[ix]->bin(iy+1), 3, ix+1, iy+1);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      // get the thrust
      const vector<double> tbins{0.7,0.8,0.85,0.9,0.95,1.0};
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // find the thrust bin
      unsigned int ithrust=0;
      for (; ithrust<tbins.size(); ++ithrust) {
        if (thrust.thrust() <= tbins[ithrust]) break;
      }
      Vector3 a1 = thrust.thrustMajorAxis();
      Vector3 a2 = thrust.thrustMinorAxis();
      // loop over the charged hadrons
      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& charged : fs.particles(Cuts::abspid==211 or Cuts::abspid==321 or Cuts::abspid==2212)) {
        double xE = charged.momentum().t()/meanBeamMom;
        double pT = sqrt(sqr(a1.dot(charged.momentum().p3()))+sqr(a2.dot(charged.momentum().p3())));
        if (charged.abspid()==211) {
          _pion[ithrust]->fill(xE, pT);
        }
        else if (charged.abspid()==321) {
          _kaon[ithrust]->fill(xE, pT);
        }
        else if (charged.abspid()==2212) {
          _proton[ithrust]->fill(xE, pT);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // scale factor including bin width in x_E
      const double fact = crossSection()/femtobarn/sumOfWeights()/0.05;
      scale(_pion, fact);
      scale(_kaon, fact);
      scale(_proton, fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _pion[6], _kaon[6], _proton[6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1718551);

}
