// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief prompt J/psi and psi(2S) at 7 TeV
  class CMS_2015_I1345023 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2015_I1345023);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // J/psi/psi(2s) histos
      for (unsigned int ipsi=0;ipsi<2;++ipsi) {
        book(_h_psi[ipsi], {0.,0.3,0.6,0.9,1.2});
        for (unsigned int iy=1;iy<5;++iy) {
          Histo1DPtr tmp;
          book(_h_psi[ipsi]->bin(iy),ipsi+1,1,iy);
        }
        book(_h_pT[ipsi][0],3+ipsi,1,1);
        book(_h_pT[ipsi][1],"TMP/psi_"+toString(ipsi),refData(5, 1, 1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for (const Particle& p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443)) {
        // prompt
        if (p.fromBottom()) continue;
        // cuts on rapidity
        const double y = p.absrap();
        if (y>1.2) continue;
        const double pT = p.perp();
        unsigned int ipsi=p.pid()/100000;
        _h_psi[ipsi]->fill(y,pT);
        _h_pT[ipsi][0]->fill(pT);
        _h_pT[ipsi][1]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 due folded rapidity
      const double factor = 0.5*crossSection() / picobarn/ sumOfWeights();
      // branching ratios to muons
      const vector<double> brs={0.05961,0.00793};
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_psi[ix], factor*brs[ix]);
        divByGroupWidth(_h_psi[ix]);
        scale(_h_pT[ix][0],brs[ix]*factor/1.2);
        scale(_h_pT[ix][1],brs[ix]*factor/1.2);
      }
      // ratio plots
      Estimate1DPtr tmp;
      book(tmp,5,1,1);
      divide(_h_pT[1][1],_h_pT[0][1],tmp);
      // convert to %age
      tmp->scale(100.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_psi[2];
    Histo1DPtr _h_pT[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2015_I1345023);

}
