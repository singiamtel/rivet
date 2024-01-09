// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  class ATLAS_2015_I1351916 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2015_I1351916);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      else if ( getOption("LMODE") == "MU" ) _mode = 2;

      LeptonFinder elecs(0.1, Cuts::pT > 25*GeV && Cuts::abseta < 4.9 && Cuts::abspid == PID::ELECTRON);
      declare(elecs, "elecs");

      LeptonFinder muons(0.1, Cuts::pT > 20*GeV && Cuts::abseta < 2.47 && Cuts::abspid == PID::MUON);
      declare(muons, "muons");

      // Book dummy histograms for heterogeneous merging
      if (_mode == 0 || _mode == 1) {
        const Estimate1D& ref_el = refData(2, 1, 2);
        book(_h["el_NCC_pos"], "_el_ncc_pos", ref_el.xEdges());
        book(_h["el_NCC_neg"], "_el_ncc_neg", ref_el.xEdges());
        book(_s["el_CC"], 2, 1, 2);

        // electron-channel only
        const Estimate1D& ref_cf = refData(3, 1, 2);
        book(_h["el_NCF_pos"], "_el_ncf_pos", ref_cf.xEdges());
        book(_h["el_NCF_neg"], "_el_ncf_neg", ref_cf.xEdges());
        book(_s["el_CF"], 3, 1, 2);
      }
      if (_mode == 0 || _mode == 2) {
        const Estimate1D& ref_mu = refData(4, 1, 2);
        book(_h["mu_NCC_pos"], "_mu_ncc_pos", ref_mu.xEdges());
        book(_h["mu_NCC_neg"], "_mu_ncc_neg", ref_mu.xEdges());
        book(_s["mu_CC"], 4, 1, 2);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get and cut on dressed leptons
      const DressedLeptons& elecs = apply<LeptonFinder>(event, "elecs").dressedLeptons();
      const DressedLeptons& muons = apply<LeptonFinder>(event, "muons").dressedLeptons();

      if (_mode == 0 || _mode == 1)  fillHistos(elecs, true);
      if (_mode == 0 || _mode == 2)  fillHistos(muons, false);
    }

    void fillHistos(const DressedLeptons& leptons, const bool doCF) {

      if (leptons.size() != 2)  return; // require exactly two leptons
      if (leptons[0].charge3() * leptons[1].charge3() > 0) return; // require opposite charge

      // Identify lepton vs antilepton
      const Particle& lpos = leptons[(leptons[0].charge3() > 0) ? 0 : 1];
      const Particle& lneg = leptons[(leptons[0].charge3() < 0) ? 0 : 1];

      string label;
      if (doCF) { // electron channel
        label = "el_N";
        const double eta1 = lpos.abseta();
        const double eta2 = lneg.abseta();
        if ( (eta1 < 2.47 && inRange(eta2, 2.5, 4.9)) || (eta2 < 2.47 && inRange(eta1, 2.5, 4.9)) ) {
          label += "CF"; // central-forward
        }
        else if (eta1 < 2.47 && eta2 < 2.47) {
          label += "CC"; // central-central
        }
        else  return; // ain't no forward-forward
      }
      else { // muon channel
        label = "mu_NCC"; // only central-central for muons
     }

      const double cosThetaStar = cosCollinsSoper(lneg, lpos);
      const double mll = (lpos.mom() + lneg.mom()).mass();
      label += cosThetaStar < 0.0?  "_neg" : "_pos";
      _h[label]->fill(mll/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSectionPerEvent() / picobarn;
      scale(_h, sf);
      if (_mode == 0 || _mode == 1) {
        asymm(_h["el_NCC_pos"], _h["el_NCC_neg"], _s["el_CC"]);
        asymm(_h["el_NCF_pos"], _h["el_NCF_neg"], _s["el_CF"]);
      }
      if (_mode == 0 || _mode == 2) {
        asymm(_h["mu_NCC_pos"], _h["mu_NCC_neg"], _s["mu_CC"]);
      }
    }


    // Cosine of the decay angle in the Collins-Soper frame
    double cosCollinsSoper(const FourMomentum& l1, const FourMomentum& l2) {
      const FourMomentum ll = l1 + l2;
      const double nom  = (l1.E() + l1.pz()) * (l2.E() - l2.pz()) - (l1.E() - l1.pz()) * (l2.E() + l2.pz());
      const double denom = ll.mass() * sqrt( sqr(ll.mass()) + sqr(ll.pt()) );
      return sign(ll.pz()) * safediv(nom, denom); // protect against division by zero, you never know...
    }

    /// @}


  private:

    /// Histograms
    map<string, Histo1DPtr> _h;
    /// Asymmetries
    map<string, Estimate1DPtr> _s;

    // option to steer the lepton channel
    size_t _mode;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2015_I1351916);

}
