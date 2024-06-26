// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {


  /// @brief onium production at 5.02 TeV
  class ATLAS_2018_I1622737 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1622737);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(Cuts::pid==443 || Cuts::pid==100443 ||
                                Cuts::pid==553 || Cuts::pid==100553 ||
                                Cuts::pid==200553), "UFS");

      // binning in y
      const vector<double> ybins={0.,0.75,1.5,2.0};
      // book histos
      for (unsigned int ix=0; ix<7; ++ix) {
        book(_h_Onium[ix], ybins);
        for (auto& b : _h_Onium[ix]->bins()) {
          book(b, ix+1, 1, b.index());
        }
        if (ix == 0 || ix == 2) {
          _axes[ix] = YODA::Axis<double>({ 8.,9.,10.,11.,12.,13.,14.,16.,18.,22.,30.,40.});
        }
        else if (ix == 1 || ix == 3) {
          _axes[ix] = YODA::Axis<double>({ 8.,10.,12.,16.,22.,40.});
        }
        else if (ix == 4) {
          _axes[ix] = YODA::Axis<double>({ 0.,2.,4.,6.,8.,10.,14.,20.,40.});
        }
        else {
          _axes[ix] = YODA::Axis<double>({ 0.,6.,10.,15.,40.});
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges[0].empty()) {
        for (unsigned int ix=0; ix<7; ++ix) {
          _edges[ix] = _h_Onium[ix]->bin(1)->xEdges();
        }
      }

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles()) {
        const double absrap = p.absrap();
        const double xp = p.perp();
        if (p.pid()==443 || p.pid()==100443) {
          // prompt/non-prompt
          const bool prompt = !p.fromBottom();
          if (p.pid()==443) {
            _h_Onium[2*prompt]->fill(absrap, disc(xp, 2*prompt));
          }
          else {
            _h_Onium[1+2*prompt]->fill(absrap, disc(xp, 1+2*prompt));
          }
        }
        else if(p.pid()==553)     _h_Onium[4]->fill(absrap,disc(xp,4));
        else if(p.pid()==100553)  _h_Onium[5]->fill(absrap,disc(xp,5));
        else if(p.pid()==200553)  _h_Onium[6]->fill(absrap,disc(xp,6));
      }
    }

    string disc(const double value, size_t ix) const {
      size_t idx = _axes[ix].index(value);
      if (0 < idx && idx <= _axes[ix].numBins()) return _edges[ix][idx-1];
      return "OTHER"s;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 due rapidity folding +/-
      const double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      const vector<double> br = {0.05961,0.00793,0.05961,0.00793,0.0248,0.0193,0.0218};
      for (unsigned int ix=0; ix<7; ++ix) {
        scale(_h_Onium[ix], factor*br[ix]);
        for(auto & hist : _h_Onium[ix]->bins()) {
          for(auto & b : hist->bins()) {
            b.scaleW(1./_axes[ix].width(b.index()));
          }
        }
        divByGroupWidth(_h_Onium[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    HistoGroupPtr<double,string> _h_Onium[7];
    YODA::Axis<double> _axes[7];
    vector<string> _edges[7];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1622737);

}
