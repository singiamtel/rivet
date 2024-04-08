// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief kaon spectra at LEP1
  class DELPHI_1995_I377487 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1995_I377487);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_K0_x , 8, 1, 1);
      book(_h_K0_xi, 9, 1, 1);
      book(_h_Ks_x ,10, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty())  _edges = _h_K0_x->xEdges();

      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==130 or Cuts::pid==310 or Cuts::abspid==323)) {
        const double xp = p.p3().mod()/meanBeamMom;
        if (abs(p.pid())==323) {
          _h_Ks_x->fill(xp);
        }
        else {
          _h_K0_x->fill(map2string(xp));
          _h_K0_xi->fill(-log(xp));
        }
      }
    }

    string map2string(const double val) const {
      const size_t idx = _axis.index(val);
      if (idx || idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_K0_x,  1./sumOfWeights());
      for(auto & b: _h_K0_x->bins()) {
        const size_t idx = b.index();
        b.scaleW(1./_axis.width(idx));
      }
      scale(_h_K0_xi, 1./sumOfWeights());
      scale(_h_Ks_x,  1./sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_K0_x;
    Histo1DPtr _h_K0_xi,_h_Ks_x;
    vector<string> _edges;
    YODA::Axis<double> _axis{ 0.00030213628351706805, 0.00048341905738733416, 0.00048343557923922857,
                              0.00026205300091453314, 0.0007702529509799692, 0.0009481642553664134,
                              0.0017226600969315643, 0.002004423179522298, 0.0016843611112658217,
                              0.0026292281438344, 0.00267627755270744, 0.0036267300396739185,
                              0.004237796021633787, 0.005212931632136056, 0.007189937374782032,
                              0.007726421785666127, 0.010282046710587495, 0.012196841637666128,
                              0.014664716763387292, 0.018701111778413465, 0.02210348200534462,
                              0.027403036058393532, 0.03380578808779788, 0.04112055882855764,
                              0.049671035882778436, 0.22518836390597363, 1.0 };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_1995_I377487);


}
