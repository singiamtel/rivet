// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Charged-particle pT and pseudorapidity spectra from pp collisions at 900 and 2360 GeV
  class CMS_2010_I845323 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2010_I845323);


    /// @{

    void init() {
      ChargedFinalState cfs(Cuts::abseta < 2.5);
      declare(cfs, "CFS");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        bool offset(en == "2360"s);

        for (int d=0; d<3; ++d) {
          for (int y=0; y<4; ++y) {
            size_t bin = 4*d+y;
            book(_h[en+"dNch_dpT"+toString(bin)], d+(offset? 4 : 1), 1, y+1);
          }
        }
        book(_h[en+"dNch_dpT_all"], 7, 1, 1+offset);
        book(_h[en+"dNch_dEta"],    8, 1, 1+offset);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    void analyze(const Event& event) {

      //charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");

      for (const Particle& p : charged.particles()) {
        //selecting only charged hadrons
        if (! PID::isHadron(p.pid())) continue;

        const double pT = p.pT();
        const double eta = p.eta();

        // The data is actually a duplicated folded distribution. This should mimic it.
        _h[_sqs+"dNch_dEta"]->fill(eta, 0.5);
        _h[_sqs+"dNch_dEta"]->fill(-eta, 0.5);
        if (fabs(eta) < 2.4 && pT > 0.1*GeV) {
          if (pT < 4.0*GeV) {
            _h[_sqs+"dNch_dpT_all"]->fill(pT/GeV, 1.0/(pT/GeV));
            if (pT < 2.0*GeV) {
              const string suff = toString(int(fabs(eta)/0.2));
              _h[_sqs+"dNch_dpT"+suff]->fill(pT/GeV);
            }
          }
        }
      }
    }


    void finalize() {
      const double normfac = 1.0/sumOfWeights(); // Normalizing to unit eta is automatic
      // The pT distributions in bins of eta must be normalized to unit eta.  This is a factor of 2
      // for the |eta| times 0.2 (eta range).
      // The pT distributions over all eta are normalized to unit eta (2.0*2.4) and by 1/2*pi*pT.
      // The 1/pT part is taken care of in the filling.  The 1/2pi is taken care of here.
      const double normpT = normfac/(2.0*0.2);
      const double normpTall = normfac/(2.0*M_PI*2.0*2.4);

      for (auto& item : _h) {
        if (item.first.find("_dEta") != string::npos) {
          scale(item.second, normfac);
        }
        else if (item.first.find("_all") != string::npos) {
          scale(item.second, normpTall);
        }
        else {
          scale(item.second, normpT);
        }
      }
    }

    /// @}


  private:

    /// @{
    map<string, Histo1DPtr> _h;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2010_I845323, CMS_2010_S8547297);

}
