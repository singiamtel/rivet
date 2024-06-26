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

      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
        iloc = 1;
      }
      else if (isCompatibleWithSqrtS(8000)) {
        iloc = 2;
      }
      else {
        throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      }
      // binning in y
      const vector<double> yedges = {0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.};
      // book histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_JPsi[ix], yedges);
        if(ix<2) book(_h_JPsi[ix+3], yedges);
        book(_h_psi2S[ix], yedges);
        for (size_t iy=1; iy < yedges.size(); ++iy) {
          if (ix == 2) {
            // total no for ratios etc
            book(_h_JPsi[ix]->bin(iy),  "TMP/JPsi_2_"+toString(iy),  refData(  iloc,1,iy));
            book(_h_psi2S[ix]->bin(iy), "TMP/psi2S_2_"+toString(iy), refData(4+iloc,1,iy));
            continue;
          }
          // prompt and non-prompt Jpsi
          book(_h_JPsi[ix]->bin(iy), iloc+2*ix,1,iy);
          // prompt and non-prompt psi(2S)
          book(_h_psi2S[ix]->bin(iy), 4+iloc+2*ix,1,iy);
          // extra Jpsi for ratios with psi2s
          if(ix<2) book(_h_JPsi[ix+3]->bin(iy), "TMP/JPsi_"+toString(ix+3)+"_"+toString(iy), refData(4+iloc+2*ix,1,iy));
        }
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
          _h_JPsi[nonPrompt  ]->fill(absrap,xp);
          _h_JPsi[2          ]->fill(absrap,xp);
          _h_JPsi[nonPrompt+3]->fill(absrap,xp);
        }
        else {
          _h_psi2S[nonPrompt]->fill(absrap,xp);
          _h_psi2S[2        ]->fill(absrap,xp);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      const double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      const vector<double> br = {0.05961,0.00793};
      // scale histos
      for (unsigned int ix=0; ix<5; ++ix) {
        scale(_h_JPsi [ix], factor*br[0]);
        divByGroupWidth(_h_JPsi[ix]);
        if(ix>2) continue;
        scale(_h_psi2S[ix], factor*br[1]);
        divByGroupWidth(_h_psi2S[ix]);
      }
      // ratios, first find CMS energy
      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
        iloc = 1;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
        iloc = 2;
      }
      else {
        throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      }

      for (unsigned int iy=1; iy<=_h_JPsi[0]->numBins(); ++iy) {
        // non-prompt J/psi percentage
        Estimate1DPtr tmp;
        book(tmp,8+iloc,1,iy);
        efficiency(_h_JPsi[1]->bin(iy),_h_JPsi[2]->bin(iy),tmp);
        tmp->scale(100.);
        // non-prompt psi2S percentage
        book(tmp,10+iloc,1,iy);
        efficiency(_h_psi2S[1]->bin(iy),_h_psi2S[2]->bin(iy),tmp);
        tmp->scale(100.);
        // prompt psi(2s)/J/psi percentage
        book(tmp,12+iloc,1,iy);
        divide(_h_psi2S[0]->bin(iy),_h_JPsi[3]->bin(iy),tmp);
        tmp->scale(100.);
        // non-prompt psi(2s)/J/psi percentage
        book(tmp,14+iloc,1,iy);
        divide(_h_psi2S[1]->bin(iy),_h_JPsi[4]->bin(iy),tmp);
        tmp->scale(100.);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_JPsi[5],_h_psi2S[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1409298);

}
