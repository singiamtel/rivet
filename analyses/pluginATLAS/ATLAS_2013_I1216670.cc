// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MPI sensitive di-jet balance variables for W->ejj or W->mujj events.
  class ATLAS_2013_I1216670 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1216670);


    /// @name Analysis methods
    /// @{

    /// Book histograms, set up projections for W and jets
    void init() {

      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      declare("MET", MissingMomentum());
      PromptFinalState ef(cuts && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");
      PromptFinalState mf(cuts && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      VetoedFinalState jet_fs;
      jet_fs.addVetoOnThisFinalState(ef);
      jet_fs.addVetoOnThisFinalState(mf);
      FastJets jets(jet_fs, JetAlg::ANTIKT, 0.4);
      declare(jets, "Jets");

      book(_h_delta_jets_n ,1, 1, 1);
      book(_h_delta_jets   ,2, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // W reco, starting with MET
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& es = apply<PromptFinalState>(event, "Elecs").particles();
      const Particles es_mtfilt = select(es, [&](const Particle& e){ return mT(e, pmiss) > 40*GeV; });
      const int iefound = closestMatchIndex(es_mtfilt, pmiss, Kin::mass, 80.4*GeV);
      const Particles& mus = apply<PromptFinalState>(event, "Muons").particles();
      const Particles mus_mtfilt = select(mus, [&](const Particle& m){ return mT(m, pmiss) > 40*GeV; });
      const int imfound = closestMatchIndex(mus_mtfilt, pmiss, Kin::mass, 80.4*GeV);

      // Event selection and ID
      if (iefound < 0 && imfound < 0) {
        MSG_DEBUG("No W's passed cuts: vetoing");
        vetoEvent;
      }
      if (iefound >= 0 && imfound >= 0) {
        MSG_DEBUG("Multiple W's passed cuts: vetoing");
        vetoEvent;
      }
      const Particle& lepton = (iefound >= 0) ? es_mtfilt[iefound] : mus_mtfilt[imfound];
      MSG_DEBUG("Event identified as W -> " << lepton.pid() << " + nu");

      // Remove jets with deltaR < 0.5 from W lepton
      const Jets all_jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20.0*GeV && Cuts::absrap < 2.8);
      const Jets jets = select(all_jets, deltaRGtr(lepton, 0.5));
      MSG_DEBUG("Overlap removal #jets = " << all_jets.size() << " -> " << jets.size());

      // Exactly two jets required
      if (jets.size() != 2)  vetoEvent;

      // Calculate analysis quantities from the two jets
      double delta_jets = (jets[0].momentum() + jets[1].momentum()).pT();
      double total_pt = jets[0].pT() + jets[1].pT();
      double delta_jets_n = delta_jets / total_pt;

      _h_delta_jets->fill(delta_jets); // Jet pT balance
      _h_delta_jets_n->fill(delta_jets_n); // Jet pT balance, normalised by scalar dijet pT
    }


    /// Finalize
    void finalize() {
      normalize(_h_delta_jets_n, 0.03);
      normalize(_h_delta_jets  , 3.0 );
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_delta_jets_n;
    Histo1DPtr _h_delta_jets;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1216670);

}
