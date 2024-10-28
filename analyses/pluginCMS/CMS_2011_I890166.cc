// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief CMS strange particle spectra (Ks, Lambda, Cascade) in pp at 900 and 7000 GeV
  ///
  /// @author Kevin Stenson
  class CMS_2011_I890166 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I890166);


    void init() {

      UnstableParticles ufs(Cuts::absrap < 2);
      declare(ufs, "UFS");

      // Particle distributions versus rapidity and transverse momentum
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        bool offset(en == "7000"s);

        book(_h[en+"dNKshort_dy"],  1, 1, 1+offset);
        book(_h[en+"dNKshort_dpT"], 2, 1, 1+offset);
        book(_h[en+"dNLambda_dy"],  3, 1, 1+offset);
        book(_h[en+"dNLambda_dpT"], 4, 1, 1+offset);
        book(_h[en+"dNXi_dy"],      5, 1, 1+offset);
        book(_h[en+"dNXi_dpT"],     6, 1, 1+offset);
        //
        book(_e[en+"LampT_KpT"],   7, 1, 1+offset);
        book(_e[en+"XipT_LampT"],  8, 1, 1+offset);
        book(_e[en+"Lamy_Ky"],     9, 1, 1+offset);
        book(_e[en+"Xiy_Lamy"],   10, 1, 1+offset);

      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    void analyze(const Event& event) {

      const UnstableParticles& parts = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : parts.particles()) {
        switch (p.abspid()) {
        case PID::K0S:
          _h[_sqs+"dNKshort_dy"]->fill(p.absrap());
          _h[_sqs+"dNKshort_dpT"]->fill(p.pT()/GeV);
          break;

        case PID::LAMBDA:
          // Lambda should not have Cascade or Omega ancestors since they should not decay. But just in case...
          if ( !( p.hasAncestorWith(Cuts::pid ==  3322) ||
                  p.hasAncestorWith(Cuts::pid == -3322) ||
                  p.hasAncestorWith(Cuts::pid ==  3312) ||
                  p.hasAncestorWith(Cuts::pid == -3312) ||
                  p.hasAncestorWith(Cuts::pid ==  3334) ||
                  p.hasAncestorWith(Cuts::pid == -3334) ) ) {
            _h[_sqs+"dNLambda_dy"]->fill(p.absrap());
            _h[_sqs+"dNLambda_dpT"]->fill(p.pT()/GeV);
          }
          break;

        case PID::XIMINUS:
          // Cascade should not have Omega ancestors since it should not decay.  But just in case...
          if ( !( p.hasAncestorWith(Cuts::pid ==  3334) ||
                  p.hasAncestorWith(Cuts::pid == -3334) ) ) {
            _h[_sqs+"dNXi_dy"]->fill(p.absrap());
            _h[_sqs+"dNXi_dpT"]->fill(p.pT()/GeV);
          }
          break;
        }
      }
    }


    void finalize() {

      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));

        divide(_h[en+"dNLambda_dpT"], _h[en+"dNKshort_dpT"], _e[en+"LampT_KpT"]);

        YODA::Histo1D denom = _h[en+"dNLambda_dpT"]->clone();
        denom.rebinXTo(_h[en+"dNXi_dpT"]->xEdges());
        divide(*_h[en+"dNXi_dpT"], denom, _e[en+"XipT_LampT"]);

        divide(_h[en+"dNLambda_dy"], _h[en+"dNKshort_dy"], _e[en+"Lamy_Ky"]);
        divide(_h[en+"dNXi_dy"], _h[en+"dNLambda_dy"], _e[en+"Xiy_Lamy"]);
      }

      const double normpT = 1.0/sumOfWeights();
      const double normy = 0.5*normpT; // Accounts for using |y| instead of y
      for (auto& item : _h) {
        if (item.first.find("_dy") != string::npos) {
          scale(item.second, normy);
        }
        else {
          scale(item.second, normpT);
        }
      }
    }


  private:

    /// @name Particle distributions versus rapidity and transverse momentum
    /// @{
    map<string,Histo1DPtr> _h;
    map<string,Estimate1DPtr> _e;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_I890166, CMS_2011_S8978280);

}
