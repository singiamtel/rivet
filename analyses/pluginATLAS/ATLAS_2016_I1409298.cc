// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {


  /// @brief ATLAS J/psi psi2s at 7 and 8 TeV
  class ATLAS_2016_I1409298 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1409298);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(UnstableParticles(Cuts::pid==443 or Cuts::pid==100443), "UFS");

      // binning in y
      const vector<double> yedges = {0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.};
      // book histos
      size_t ie = 0; _ih = 10;
      for (double eVal : allowedEnergies()) {

        if (isCompatibleWithSqrtS(eVal))  _ih = ie;

        for (unsigned int ix=0; ix<3; ++ix) {
          book(_h_JPsi[ie][ix], yedges);
          if(ix<2) book(_h_JPsi[ie][ix+3], yedges);
          book(_h_psi2S[ie][ix], yedges);
          for (size_t iy=1; iy < yedges.size(); ++iy) {
            if (ix == 2) {
              // total no for ratios etc
              book(_h_JPsi[ie][ix]->bin(iy),
                   "TMP/JPsi_2_"+toString(iy)+"_"+toString(ie),
                   refData(1+ie,1,iy));
              book(_h_psi2S[ie][ix]->bin(iy),
                   "TMP/psi2S_2_"+toString(iy)+"_"+toString(ie),
                   refData(5+ie,1,iy));
              continue;
            }
            // prompt and non-prompt Jpsi
            book(_h_JPsi[ie][ix]->bin(iy), 1+ie+2*ix,1,iy);
            // prompt and non-prompt psi(2S)
            book(_h_psi2S[ie][ix]->bin(iy), 5+ie+2*ix,1,iy);
            // extra Jpsi for ratios with psi2s
            if(ix<2) {
              book(_h_JPsi[ie][ix+3]->bin(iy),
                   "TMP/JPsi_"+toString(ix+3)+"_"+toString(iy)+"_"+toString(ie),
                   refData(5+ie+2*ix,1,iy));
            }
          }
        }
        ++ie;
      }
      if (_ih == 10 && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles()) {
        // prompt/non-prompt
        bool nonPrompt = p.fromBottom();
        const double absrap = p.absrap();
        const double xp = p.perp();
        if (p.pid()==443) {
          _h_JPsi[_ih][nonPrompt  ]->fill(absrap,xp);
          _h_JPsi[_ih][2          ]->fill(absrap,xp);
          _h_JPsi[_ih][nonPrompt+3]->fill(absrap,xp);
        }
        else {
          _h_psi2S[_ih][nonPrompt]->fill(absrap,xp);
          _h_psi2S[_ih][2        ]->fill(absrap,xp);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      const double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      const vector<double> br = {0.05961,0.00793};
      for (size_t ie=0; ie<2; ++ie) {
        // scale histos
        for (size_t ix=0; ix<5; ++ix) {
          scale(_h_JPsi[ie][ix], factor*br[0]);
          divByGroupWidth(_h_JPsi[ie][ix]);
          if (ix>2) continue;
          scale(_h_psi2S[ie][ix], factor*br[1]);
          divByGroupWidth(_h_psi2S[ie][ix]);
        }
        // ratios
        for (size_t iy=1; iy<=_h_JPsi[ie][0]->numBins(); ++iy) {
          // non-prompt J/psi percentage
          Estimate1DPtr tmp;
          book(tmp,9+ie,1,iy);
          efficiency(_h_JPsi[ie][1]->bin(iy),_h_JPsi[ie][2]->bin(iy),tmp);
          tmp->scale(100.);
          // non-prompt psi2S percentage
          book(tmp,11+ie,1,iy);
          efficiency(_h_psi2S[ie][1]->bin(iy),_h_psi2S[ie][2]->bin(iy),tmp);
          tmp->scale(100.);
          // prompt psi(2s)/J/psi percentage
          book(tmp,13+ie,1,iy);
          divide(_h_psi2S[ie][0]->bin(iy),_h_JPsi[ie][3]->bin(iy),tmp);
          tmp->scale(100.);
          // non-prompt psi(2s)/J/psi percentage
          book(tmp,15+ie,1,iy);
          divide(_h_psi2S[ie][1]->bin(iy),_h_JPsi[ie][4]->bin(iy),tmp);
          tmp->scale(100.);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_JPsi[2][5], _h_psi2S[2][3];
    size_t _ih;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1409298);

}
