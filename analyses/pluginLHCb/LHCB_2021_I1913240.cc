// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {

  /// @brief nuclear modification factor and prompt charged particle production in pPb and pp at $\sqrt{s_{NN}} = 5$ TeV
  class LHCB_2021_I1913240 : public Analysis {

  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2021_I1913240);

    void init() {
      // Check collision system of input file
      const ParticlePair &beam_pair = beams();
      if (beam_pair.first.pid() == PID::PROTON &&
          beam_pair.second.pid() == PID::PROTON) {
        _beamConfig = "pp";
      }
      else if (beam_pair.first.pid() == PID::PROTON &&
          beam_pair.second.pid() == PID::LEAD) {
        _beamConfig = "pPb";
      }
      if (_beamConfig == "" && !merging()) {
        throw BeamError("Invalid beam configuration for " + name() + "\n");
      }

      // Link histogram to corresponding table of YODA file from HEPData
      book(_charged_part_eta_pt_hist[0], {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 4.8},
        {"d01-x01-y01", "d01-x01-y02", "d01-x01-y03", "d01-x01-y04", "d01-x01-y05", "d01-x01-y06"});
      // forward and backward histos for pPb below
      book(_charged_part_eta_pt_hist[1], {1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.3},
           {"d02-x01-y01", "d02-x01-y02", "d02-x01-y03", "d02-x01-y04", "d02-x01-y05", "d02-x01-y06"});
      book(_charged_part_eta_pt_hist[2], {-5.2, -4.8, -4.5, -4.0, -3.5, -3.0, -2.5},
           {"d03-x01-y01", "d03-x01-y02", "d03-x01-y03", "d03-x01-y04", "d03-x01-y05", "d03-x01-y06"});

      // Select primary (i.e. prompt long-lived) charged particles
      declare(ALICE::PrimaryParticles(
                  (Cuts::absetaIn(1.6, 5.2)) &&
                  (Cuts::abscharge == 1)),
              "PPs");
    }

    void analyze(const Event &ev) {
      // Apply PrimaryParticles projection
      const Particles prompt_charged_parts = apply<ALICE::PrimaryParticles>(ev, "PPs").particles();
      for (const Particle &part : prompt_charged_parts) {
        if (part.momentum().p() < 2 * GeV)  continue;

        // Double sample size in symmetric proton-proton case by taking absolute value
        // and then scale by 0.5 due to LHCb being one-sided
        double eta;
        if (_beamConfig == "pp")  eta = part.abseta();
        else                   eta = part.eta();
        double pt = part.pT();

        if (_beamConfig == "pp") {
          _charged_part_eta_pt_hist[0]->fill(eta, pt, 0.5);
        }
        else {
          _charged_part_eta_pt_hist[1]->fill(eta, pt);
          _charged_part_eta_pt_hist[2]->fill(eta, pt);
        }
      }
    }

    void finalize() {
      // Compute scale factor with inelastic cross-section from input file and sum of
      // weights (corresponds to number of events in input file)
      divByGroupWidth(_charged_part_eta_pt_hist);
      scale(_charged_part_eta_pt_hist, crossSection() / millibarn / sumOfWeights());

      // if both histo groups have filled bins then compute R_pPb
      Estimate1DPtr _ratio;
      for (size_t i = 1; i <= 5; ++i) {
        book(_ratio, 4, 1, i);
        YODA::Histo1D numer = _charged_part_eta_pt_hist[1]->bin(i+1)->clone();
        YODA::Histo1D denom = _charged_part_eta_pt_hist[0]->bin(i)->clone();
        numer.rebinXTo(_ratio->xEdges());
        denom.rebinXTo(_ratio->xEdges());
        divide(numer, denom, _ratio);
        _ratio->scale(1./208.);
        book(_ratio, 5, 1, i);
        numer = _charged_part_eta_pt_hist[2]->bin(i+1)->clone();
        denom = _charged_part_eta_pt_hist[0]->bin(7-i)->clone();
        numer.rebinXTo(_ratio->xEdges());
        denom.rebinXTo(_ratio->xEdges());
        divide(numer, denom, _ratio);
        _ratio->scale(1./208.);
      }

    }

    string _beamConfig("");
    Histo1DGroupPtr _charged_part_eta_pt_hist[3];
  };

  RIVET_DECLARE_PLUGIN(LHCB_2021_I1913240);
}
