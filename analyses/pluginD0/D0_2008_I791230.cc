// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief D0 Run II measurement of W charge asymmetry
  ///
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_I791230 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2008_I791230);


    /// @name Analysis methods
    /// @{

    // Book histograms and set up projections
    void init() {

      // Projections
      declare("MET", MissingMomentum());
      LeptonFinder ef(0.2, Cuts::abseta < 5 && Cuts::pT > 25*GeV && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");

      // Histograms (temporary +- charge histos and scatters to store the calculated asymmetries)
      for (size_t pmindex = 0; pmindex < 2; ++pmindex) {
        const string suffix = (pmindex == 0) ? "plus" : "minus";
        book(_hs_dsigpm_deta_25_35[pmindex] ,"TMP/dsigpm_deta_25_35_" + suffix, refData(1, 1, 1));
        book(_hs_dsigpm_deta_35[pmindex] ,"TMP/dsigpm_deta_35_" + suffix, refData(1, 1, 2));
        book(_hs_dsigpm_deta_25[pmindex] ,"TMP/dsigpm_deta_25_" + suffix, refData(1, 1, 3));
      }
      book(_h_asym1, 1, 1, 1);
      book(_h_asym2, 1, 1, 2);
      book(_h_asym3, 1, 1, 3);
    }


    /// Do the analysis
    void analyze(const Event & event) {

      // W reco, starting with MET
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int ifound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      if (ifound < 0) {
        MSG_DEBUG("No W candidates found: vetoing");
        vetoEvent;
      }

      // Get the e+- momentum, and an effective charge including the eta sign
      const Particle& e = es[ifound];
      const int chg_e = sign(e.eta()) * sign(e.charge());
      assert(chg_e == 1 || chg_e == -1);
      MSG_TRACE("Charged lepton sign = " << chg_e);

      // Fill histos with appropriate +- indexing
      const size_t pmindex = (chg_e > 0) ? 0 : 1;
      if (e.Et() < 35*GeV) _hs_dsigpm_deta_25_35[pmindex]->fill(e.abseta());
      else _hs_dsigpm_deta_35[pmindex]->fill(e.abseta());
      _hs_dsigpm_deta_25[pmindex]->fill(e.abseta());
    }


    /// @brief Finalize
    ///
    /// Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each ET region
    void finalize() {
      asymm(_hs_dsigpm_deta_25_35[0], _hs_dsigpm_deta_25_35[1], _h_asym1);
      asymm(_hs_dsigpm_deta_35[0],    _hs_dsigpm_deta_35[1],    _h_asym2);
      asymm(_hs_dsigpm_deta_25[0],    _hs_dsigpm_deta_25[1],    _h_asym3);
      _h_asym1->scale(100.);
      _h_asym2->scale(100.);
      _h_asym3->scale(100.);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _hs_dsigpm_deta_25_35[2], _hs_dsigpm_deta_35[2], _hs_dsigpm_deta_25[2];
    Estimate1DPtr _h_asym1, _h_asym2, _h_asym3;
    /// @}

  };


  RIVET_DECLARE_ALIASED_PLUGIN(D0_2008_I791230, D0_2008_S7837160);

}
