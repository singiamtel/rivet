// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > c cbar at 10.52 GeV
  class CLEO_1988_I23509 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1988_I23509);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==411  ||
                                Cuts::abspid==421  ||
                                Cuts::abspid==431  ||
                                Cuts::abspid==4122 ||
                                Cuts::abspid==4232 ||
                                Cuts::abspid==4132 ||
                                Cuts::abspid==4332), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],2+ix,1,1);
      }
      book(_h[2],"TMP/h_mu",refData<YODA::BinnedEstimate<string>>(2,1,2));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // mu+mu- + photons
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        _h[2]->fill("10.5"s);
      }
      else { // everything else
        const FinalState& ufs = apply<FinalState>(event, "UFS");
        if(!ufs.particles().empty()) {
          _h[0]->fill("10.5"s);
          _h[1]->fill("10.5"s);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection()/sumOfWeights()/nanobarn);
      BinnedEstimatePtr<string> tmp;
      book(tmp, 2, 1, 2);
      divide(_h[0], _h[2],tmp);
      book(tmp, 3, 1, 2);
      divide(_h[1], _h[2],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1988_I23509);

}
