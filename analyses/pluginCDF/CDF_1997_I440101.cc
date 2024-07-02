// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief J/psi and psi(2S) at 1.8 TeV
  class CDF_1997_I440101 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1997_I440101);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_total[ix], 1, 1, 1+ix);
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_h_psi[ix][iy], 2+2*ix+iy, 1, 1);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        _edges.resize(2);
        for(unsigned int ix=0;ix<2;++ix) {
          _edges[ix].resize(2);
          for(unsigned int iy=0;iy<2;++iy)
            _edges[ix][iy] = _h_psi[ix][iy]->xEdges();
        }
      }

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443)) {
        const double abseta = p.abseta();
        const double xp = p.perp();
          if (xp<5. || abseta>0.6) continue;
          bool prompt = !p.fromBottom();
          unsigned int ipsi=p.pid()/100000;
          _h_total[ipsi]->fill("< 0.6"s);
          _h_psi[prompt][ipsi]->fill(disc(xp, prompt, ipsi));
      }
    }

    string disc(const double value, unsigned int prompt, unsigned int ipsi) const {
      size_t idx = _axis[ipsi].index(value);
      if (0 < idx && idx <= _axis[ipsi].numBins())
        return _edges[prompt][ipsi][idx-1];
      return "OTHER"s;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalisation factor
      const double factor = crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      const vector<double> br = {0.05961,0.00793};
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_total[ix], factor*br[ix]);
        for (unsigned int iy=0; iy<2; ++iy) {
          scale(_h_psi[ix][iy], factor*br[iy]);
          for(auto & b :_h_psi[ix][iy]->bins()) {
            b.scaleW(1./_axis[iy].width(b.index()));
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_total[2], _h_psi[2][2];
    YODA::Axis<double> _axis[2] = { YODA::Axis<double>({5.,5.5,6.,6.5,7.,8.,9.,10.,12.,14.,17.,20.}),
                                    YODA::Axis<double>({5.,6.,7.,9.,12.,16.}) };
    vector<vector<vector<string>>> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CDF_1997_I440101);

}
