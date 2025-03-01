// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {

  class LHCB_2019_I1720413 : public Analysis {

  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1720413);

    void init() {

      size_t ih = 1;
      vector<double> ybins = { -4.5, -3.5, -2.5, 1.5, 2.5, 3.5 };
      for (const string& label : vector<string>{"b_plus", "b_0", "lambda_b_0"}) {
        book(_g[label], ybins);
        _g[label]->maskBin(3);
        size_t idx = 1;
        for (auto& b : _g[label]->bins()) {
          book(b, ih, 1, idx++);
        }
        ++ih;
      }

      book(_h["b_plus"], 4, 1, 1);
      book(_h["b_0"], 4, 1, 2);
      book(_h["lambda_b_0"], 4, 1, 3);

      book(_e["b_plus_diff"], 5, 1, 1);
      book(_e["b_0_diff"], 5, 1, 2);
      book(_e["lambda_b_0_diff"], 5, 1, 3);

      book(_e["b_plus"], 6, 1, 1);
      book(_e["b_0"], 6, 1, 2);
      book(_e["lambda_b_0"], 6, 1, 3);

      // Define histograms with one rapidity interval to compute
      // transverse-momentum-integrated forward-to-backward ratio
      const vector<double> forw_y_ints = {2.0, 20.0};
      book(_h["b_plus_fwd"], "_b_plus_fwd", forw_y_ints);
      book(_h["b_plus_bwd"], "_b_plus_bwd", forw_y_ints);
      book(_h["b_0_fwd"], "_b_0_fwd", forw_y_ints);
      book(_h["b_0_bwd"], "_b_0_bwd", forw_y_ints);
      book(_h["lambda_b_0_fwd"], "_lambda_b_0_fwd", forw_y_ints);
      book(_h["lambda_b_0_bwd"], "_lambda_b_0_bwd", forw_y_ints);

      // Select B+, B0 and Lambdab0 hadrons and their antiparticles
      declare(
          UnstableParticles(
              ((Cuts::abspid == PID::BPLUS) ||
               (Cuts::abspid == PID::B0) ||
               (Cuts::abspid == PID::LAMBDAB)) &&
              (Cuts::ptIn(2 * GeV, 20 * GeV)) &&
              (Cuts::rapIn(-4.5, -2.5) || Cuts::rapIn(1.5, 3.5))), "UPs");
    }

    void analyze(const Event &event) {
      // Apply UnstableParticles projection
      const UnstableParticles &unst_parts = apply<UnstableParticles>(event, "UPs");
      for (const Particle &part : unst_parts.particles()) {
        const double y = part.rapidity();
        const double pt = part.pT();

        if (part.abspid () == PID::BPLUS) {
          _g["b_plus"]->fill(y, pt/GeV);
          _h["b_plus"]->fill(y);
          if (fabs(y) < 2.5 || fabs(y) > 3.5)  continue;
          if (y>0.)  _h["b_plus_fwd"]->fill(pt/GeV);
          else       _h["b_plus_bwd"]->fill(pt/GeV);
        }
        else if (part.abspid() == PID::B0) {
          _g["b_0"]->fill(y, pt/GeV);
          _h["b_0"]->fill(y);
          if (fabs(y) < 2.5 || fabs(y) > 3.5)  continue;
          if (y>0.)  _h["b_0_fwd"]->fill(pt/GeV);
          else       _h["b_0_bwd"]->fill(pt/GeV);
        }
        else if (part.abspid() == PID::LAMBDAB) {
          _g["lambda_b_0"]->fill(y, pt/GeV);
          _h["lambda_b_0"]->fill(y);
          if (fabs(y) < 2.5 || fabs(y) > 3.5)  continue;
          if (y>0.)  _h["lambda_b_0_fwd"]->fill(pt/GeV);
          else       _h["lambda_b_0_bwd"]->fill(pt/GeV);
        }
      }
    }

    void finalize() {
      // Compute scale factor with inelastic cross-section from input file and sum of
      // weights (corresponds to number of events in input file)
      const double sf = crossSection() / microbarn / sumOfWeights();
      scale(_g, sf);
      scale(_h, sf);

      // Compute forward-to-backward ratio with correct indices from HistoGroup
      for (auto& grp : _g) {
        divide(grp.second->bin(5), grp.second->bin(2), _e[grp.first + "_diff"s]);
        divide(_h[grp.first+"_fwd"s], _h[grp.first+"_bwd"s], _e[grp.first]);
      }
    }

  private:

    map<string,Histo1DPtr> _h;
    map<string,Estimate1DPtr> _e;
    map<string,Histo1DGroupPtr> _g;
  };

  RIVET_DECLARE_PLUGIN(LHCB_2019_I1720413);
}
