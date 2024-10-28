// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Ratio of energy in -6.6 < eta < -5.2 for events with a charged-particle jet
  class CMS_2013_I1218372 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2013_I1218372);


    void init() {

      // gives the range of eta and min pT for the final state from which I get the jets
      FastJets jetpro (ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT >= 0.3*GeV), JetAlg::ANTIKT, 0.5);
      declare(jetpro, "Jets");

      // skip Neutrinos and Muons
      VetoedFinalState fsv(FinalState(Cuts::etaIn(-7.0, -4.0)));
      fsv.vetoNeutrinos();
      fsv.addVetoPairId(PID::MUON);
      declare(fsv, "fsv");

      // for the hadron level selection
      VetoedFinalState sfsv;
      sfsv.vetoNeutrinos();
      sfsv.addVetoPairId(PID::MUON);
      declare(sfsv, "sfsv");

      size_t id = 0;
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        ++id;

        //counters
        book(_c[en+"pass"], "_pass"+en);
        book(_c[en+"incl"], "_incl"+en);

        // Temporary histograms to fill the energy flow for leading jet events.
        // Ratios are calculated in finalize().
        book(_e[en+"ratio"], id, 1, 1);
        book(_h[en+"num"], "TMP/eflow_jet"+en,  refData(id, 1, 1)); // Leading jet energy flow in pt
        book(_h[en+"den"], "TMP/number_jet"+en, refData(id, 1, 1)); // Number of events in pt
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Skip if the event is empty
      const FinalState& fsv = apply<FinalState>(event, "fsv");
      if (fsv.empty()) vetoEvent;

      // ====================== Minimum Bias selection

      const FinalState& sfsv = apply<FinalState>(event, "sfsv");
      Particles parts = sfsv.particles(cmpMomByRap);
      if (parts.empty()) vetoEvent;

      // find dymax
      double dymax = 0;
      int gap_pos  = -1;
      for (size_t i = 0; i < parts.size()-1; ++i) {
        double dy = parts[i+1].rap() - parts[i].rap();
        if (dy > dymax) {
          dymax = dy;
          gap_pos = i;
        }
      }

      // calculate mx2 and my2
      FourMomentum xmom;
      for (int i=0; i<=gap_pos; ++i) {
        xmom += parts[i].mom();
      }
      double mx2 = xmom.mass2();
      if (mx2<0) vetoEvent;

      FourMomentum ymom;
      for (size_t i=gap_pos+1; i<parts.size(); ++i) {
        ymom += parts[i].mom();
      }
      double my2 = ymom.mass2();
      if (my2<0) vetoEvent;

      // calculate xix and xiy and xidd
      double xix  = mx2 / sqr(sqrtS());
      double xiy  = my2 / sqr(sqrtS());
      double xidd = mx2*my2 / sqr(sqrtS()*0.938*GeV);

      // combine the selection: xi cuts
      bool passedHadronCuts = false;
      if (_sqs=="900"s  && (xix > 0.1  || xiy > 0.4 || xidd > 0.5))  passedHadronCuts = true;
      if (_sqs=="2760"s && (xix > 0.07 || xiy > 0.2 || xidd > 0.5))  passedHadronCuts = true;
      if (_sqs=="7000"s && (xix > 0.04 || xiy > 0.1 || xidd > 0.5))  passedHadronCuts = true;
      if (!passedHadronCuts) vetoEvent;

      //  ============================== MINIMUM BIAS EVENTS

      // loop over particles to calculate the energy
      _c[_sqs+"pass"]->fill();

      for (const Particle& p : fsv.particles()) {
        if (-5.2 > p.eta() && p.eta() > -6.6) _c[_sqs+"incl"]->fill(p.E()/GeV);
      }

      //  ============================== JET EVENTS

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets& jets = jetpro.jetsByPt(Cuts::pT > 1.0*GeV);
      if (jets.size() < 1) vetoEvent;

      if (jets[0].abseta() < 2.0) {
        _h[_sqs+"den"]->fill(jets[0].pT()/GeV);

        // energy flow
        for (const Particle& p : fsv.particles()) {
          if (p.eta() > -6.6 && p.eta() < -5.2) {  // ask for the CASTOR region
            _h[_sqs+"num"]->fill(jets[0].pT()/GeV, p.E()/GeV);
          }
        }
      }

    }// analysis

    void finalize() {
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (_c[en+"incl"]->sumW()) scale(_h[en+"num"], *_c[en+"pass"] / *_c[en+"incl"]);
        divide(_h[en+"num"], _h[en+"den"], _e[en+"ratio"]);
      }
    }

  private:
    // counters
    map<string,CounterPtr> _c;
    map<string,Histo1DPtr> _h;
    map<string,Estimate1DPtr> _e;

    string _sqs = "";
  };


  RIVET_DECLARE_PLUGIN(CMS_2013_I1218372);

}
