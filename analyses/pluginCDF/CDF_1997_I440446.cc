// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi from chi_c decays
  class CDF_1997_I440446 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1997_I440446);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      const YODA::BinnedEstimate<int>& ref = refData<YODA::BinnedEstimate<int>>(1,1,1);
      book(_h_total[0], "TMP/h_chi", ref);
      book(_h_total[1], "TMP/h_psi", ref);
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_psi[ix], 2+ix, 1, 1);
      }
      _axis = YODA::Axis<double>({5.,5.5,6.,6.5,7.,8.,9.,10.,12.,14.,17.,20.});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_psi[0]->xEdges();

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // first J/psi for denominator
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
        if (p.fromBottom()) continue;
        const double abseta = p.abseta();
        const double xp = p.perp();
        if (xp<4. || abseta>0.6) continue;
        _h_total[1]->fill(1800);
        // from those for higher charmonium
        Particle parent = p.parents()[0];
        if ( parent.pid()==100443 || parent.pid()==20443 || parent.pid()==445) continue;
        _h_psi[0]->fill(disc(xp));
      }
      // chi_1 and chi_2 for numerator
      for (const Particle& p : ufs.particles(Cuts::pid==20443 || Cuts::pid==100443 || Cuts::pid==445)) {
        if(p.fromBottom()) continue;
        Particle Jpsi;
        bool found(false);
        for (const Particle & child : p.children()) {
          if (child.pid()==443) {
            found = true;
            Jpsi=child;
          }
        }
        if(!found) continue;
        double abseta=Jpsi.abseta();
        double xp = Jpsi.perp();
        if (xp<4. || abseta>0.6) continue;
        if (p.pid()==100443) {
          _h_psi[2]->fill(disc(xp));
        }
        else {
          _h_psi[1]->fill(disc(xp));
          _h_total[0]->fill(1800);
        }
      }
    }

    string disc(const double value) {
      size_t idx = _axis.index(value);
      if (0 < idx && idx <= _axis.numBins())  return _edges[idx-1];
      return "OTHER"s;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalisation factor
      const double br = 0.05961;
      const double factor = br*crossSection()/nanobarn/sumOfWeights();
      scale(_h_psi,factor);
      for(unsigned int ix=0;ix<3;++ix) {
        for(auto & b : _h_psi[ix]->bins())
          b.scaleW(1./_axis.width(b.index()));
      }
      BinnedEstimatePtr<int> tmp;
      book(tmp, 1, 1, 1);
      efficiency(_h_total[0], _h_total[1], tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _h_total[2];
    BinnedHistoPtr<string> _h_psi[3];
    YODA::Axis<double> _axis;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CDF_1997_I440446);

}
