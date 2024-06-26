// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief charm spectra at 10.55
  class CLEO_1988_I22954 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1988_I22954);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==411 ||
                                Cuts::abspid==421 ||
                                Cuts::abspid==413 ||
                                Cuts::abspid==423 ||
                                Cuts::abspid==431 ||
                                Cuts::abspid==4122), "UFS");
      // histos
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h_spect[ix],1+2*ix,1,1);
        book(_h_scale[ix],2+2*ix,1,1);
        book(_h_sigma[ix],14,1,1+ix);
        if (ix==3) continue;
        book(_h_spect[4+ix], 9+ix,1,1);
      }
      book(_h_scale[4],12,1,1);
      for (unsigned int ix=0;ix<2;++ix) {
	      for (unsigned int iy=0;iy<7;++iy) {
          book(_h_sigma_br[ix][iy],13,1+ix,1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        double Emax = 0.5*sqrtS();
        double pmax = sqrt(sqr(Emax)-p.mass2());
        double xp = (p.E()+p.p3().mod())/(Emax+pmax);
        if (p.abspid()==413) {
          _h_spect[0]->fill(xp);
          _h_scale[0]->fill(xp);
          _h_spect[1]->fill(xp);
          _h_scale[1]->fill(xp);
          _h_sigma[2]->fill("10.55"s);
          for (unsigned int ix=0; ix<2; ++ix) {
            _h_sigma_br[ix][2]->fill(xp);
          }
          for (unsigned int ix=0; ix<2; ++ix) {
            _h_sigma_br[ix][3]->fill(xp);
          }
        }
        else if (p.abspid()==421) {
          _h_spect[2]->fill(xp);
          _h_scale[2]->fill(xp);
          _h_sigma[0]->fill("10.55"s);
          for (unsigned int ix=0; ix<2; ++ix) {
            _h_sigma_br[ix][0]->fill(xp);
          }
        }
        else if (p.abspid()==411) {
          _h_spect[3]->fill(xp);
          _h_scale[3]->fill(xp);
          _h_sigma[1]->fill("10.55"s);
          for (unsigned int ix=0; ix<2; ++ix) {
            _h_sigma_br[ix][1]->fill(xp);
          }
        }
        else if (p.abspid()==431) {
          _h_spect[4]->fill(xp);
          for (unsigned int ix=0;ix<2;++ix) {
            _h_sigma_br[ix][5]->fill(xp);
          }
        }
        else if (p.abspid()==4122) {
          _h_spect[5]->fill(xp);
          for (unsigned int ix=0;ix<2;++ix) {
            _h_sigma_br[ix][6]->fill(xp);
          }
        }
        else if (p.abspid()==423) {
          _h_spect[6]->fill(xp);
          _h_scale[4]->fill(xp);
          _h_sigma[3]->fill("10.55"s);
          for (unsigned int ix=0;ix<2;++ix) {
            _h_sigma_br[ix][4]->fill(xp);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const vector<double> br = {0.677*0.03947,0.677*0.0822,0.03947,0.0938,0.045,0.0628,0.03947};
      const double fact1 = crossSection()/picobarn/sumOfWeights();
      const double fact2 = crossSection()/nanobarn/sumOfWeights()*sqr(sqrtS());
      for (unsigned int ix=0; ix<7; ++ix) {
        scale(_h_spect[ix],fact1*br[ix]);
        if (ix>4) continue;
        scale(_h_scale[ix],fact2);
        if (ix==4) continue;
        scale(_h_sigma[ix],fact1*1e-3);
      }
      const double width[2]={1.,0.5};
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_sigma_br[ix][0],width[ix]*fact1*br[2]);
        scale(_h_sigma_br[ix][1],width[ix]*fact1*br[3]);
        scale(_h_sigma_br[ix][2],width[ix]*fact1*br[0]);
        scale(_h_sigma_br[ix][3],width[ix]*fact1*br[1]);
        scale(_h_sigma_br[ix][4],width[ix]*fact1*br[6]);
        scale(_h_sigma_br[ix][5],width[ix]*fact1*br[4]);
        scale(_h_sigma_br[ix][6],width[ix]*fact1*br[5]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect[7],_h_scale[5],_h_sigma_br[2][7];
    BinnedHistoPtr<string> _h_sigma[4];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CLEO_1988_I22954);

}
