// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D meson production at 7 TeV
  class ATLAS_2016_I1408878 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1408878);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==411 or Cuts::abspid==413 or Cuts::abspid==431), "UFS");
      // histograms
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_tot[ix], "TMP/_total_"+toString(ix+1), refData(1, 1, 1+ix));
        book(_e_tot[ix], 1, 1, 1+ix);
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_pT[ix], 2, 1, 1+ix);
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_h_y[iy][ix], 3+ix, 1, 1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if (p.children().size()==1) continue;
        const double abseta = p.abseta();
        if (abseta>2.1) continue;
        const double pT = p.perp();
        if (pT<3.5 || pT>100.) continue;

        unsigned int itype=0;
        if (p.abspid()==411) itype=1;
        else if(p.abspid()==431) itype=2;
        _h_tot[itype]->fill(pT);
        if (itype<2) {
          _h_pT[itype]->fill(pT);
          if (pT<20.)  _h_y[itype][0]->fill(abseta);
          else         _h_y[itype][1]->fill(abseta);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double factor = crossSection()/microbarn/sumOfWeights();
      for (unsigned int ix=0;ix<2;++ix) {
        scale(_h_tot[ix],factor);
        scale(_h_pT[ix],factor);
        scale(_h_y[ix][0],factor);
        scale(_h_y[ix][1],1000.*factor);
      }
      scale(_h_tot[2],factor);
      for (unsigned int ix=0;ix<3;++ix) {
        barchart(_h_tot[ix], _e_tot[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2],_h_y[2][2],_h_tot[3];
    Estimate1DPtr _e_tot[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1408878);

}
