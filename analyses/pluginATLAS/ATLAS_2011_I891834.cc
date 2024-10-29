// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Calo-based underlying event at 900 GeV and 7 TeV in ATLAS
  ///
  /// @author Jinlong Zhang
  class ATLAS_2011_I891834 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I891834);


    void init() {
      const FinalState fs500(Cuts::abseta < 2.5 && Cuts::pT >=  500*MeV);
      declare(fs500, "FS500");
      const FinalState fslead(Cuts::abseta < 2.5 && Cuts::pT >=  1.0*GeV);
      declare(fslead, "FSlead");

      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        bool is7TeV(en == "7000");

        // N profiles, 500 MeV pT cut
        book(_h[en+"N_transverse"], 1+is7TeV, 1, 1);
        // pTsum profiles, 500 MeV pT cut
        book(_h[en+"ptsum_transverse"], 3+is7TeV, 1, 1);
        // N vs. Delta(phi) profiles, 500 MeV pT cut
        book(_h[en+"N_vs_dPhi_1"], 13+is7TeV, 1, 1);
        book(_h[en+"N_vs_dPhi_2"], 13+is7TeV, 1, 2);
        book(_h[en+"N_vs_dPhi_3"], 13+is7TeV, 1, 3);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    void analyze(const Event& event) {
      // Require at least one cluster in the event with pT >= 1 GeV
      const FinalState& fslead = apply<FinalState>(event, "FSlead");
      if (fslead.size() < 1)  vetoEvent;

      // These are the particles  with pT > 500 MeV
      const FinalState& chargedNeutral500 = apply<FinalState>(event, "FS500");

      // Identify leading object and its phi and pT
      Particles particles500 = chargedNeutral500.particlesByPt();
      Particle p_lead = particles500[0];
      const double philead = p_lead.phi();
      const double etalead = p_lead.eta();
      const double pTlead  = p_lead.pT();
      MSG_DEBUG("Leading object: pT = " << pTlead << ", eta = " << etalead << ", phi = " << philead);

      // Iterate over all > 500 MeV particles and count particles and scalar pTsum in the three regions
      vector<double> num500(3, 0), ptSum500(3, 0.0);
      // Temporary histos that bin N in dPhi.
      // NB. Only one of each needed since binnings are the same for the energies and pT cuts
      YODA::Histo1D htemp(_h[_sqs+"N_vs_dPhi_1"]->xEdges(),"tmp");
      for (const Particle& p : particles500) {
        const double pT = p.pT();
        const double dPhi = deltaPhi(philead, p.phi());
        const int ir = region_index(dPhi);
        num500[ir] += 1;
        ptSum500[ir] += pT;

        // Fill temp histos to bin N in dPhi
        if (p.genParticle() != p_lead.genParticle()) { // We don't want to fill all those zeros from the leading track...
          htemp.fill(dPhi);
        }
      }


      // Now fill underlying event histograms
      // The densities are calculated by dividing the UE properties by dEta*dPhi
      // -- each region has a dPhi of 2*PI/3 and dEta is two times 2.5
      const double dEtadPhi = (2*2.5 * 2*PI/3.0);
      _h[_sqs+"N_transverse"]->fill(pTlead/GeV,  num500[1]/dEtadPhi);
      _h[_sqs+"ptsum_transverse"]->fill(pTlead/GeV, ptSum500[1]/GeV/dEtadPhi);

      // Update the "proper" dphi profile histograms
      // Note that we fill dN/dEtadPhi: dEta = 2*2.5, dPhi = 2*PI/nBins
      // The values tabulated in the note are for an (undefined) signed Delta(phi) rather than
      // |Delta(phi)| and so differ by a factor of 2: we have to actually norm for angular range = 2pi
      for (const auto& b : htemp.bins()) {
        double mean = b.xMid(), value = 0.;
        if (b.numEntries() > 0) {
          mean = b.xMean();
          value = b.sumW()/b.xWidth()/10.0;
        }
        if (pTlead/GeV >= 1.0) _h[_sqs+"N_vs_dPhi_1"]->fill(mean, value);
        if (pTlead/GeV >= 2.0) _h[_sqs+"N_vs_dPhi_2"]->fill(mean, value);
        if (pTlead/GeV >= 3.0) _h[_sqs+"N_vs_dPhi_3"]->fill(mean, value);
      }

    }


    // Little helper function to identify Delta(phi) regions
    inline int region_index(double dphi) {
      assert(inRange(dphi, 0.0, PI, CLOSED, CLOSED));
      if (dphi < PI/3.0) return 0;
      if (dphi < 2*PI/3.0) return 1;
      return 2;
    }

  private:

    string _sqs = "";

    map<string,Profile1DPtr> _h;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_I891834, ATLAS_2011_S8994773);

}
