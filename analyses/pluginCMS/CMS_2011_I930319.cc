// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Forward energy flow in MB and dijet events at 0.9 and 7 TeV
  class CMS_2011_I930319 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I930319);


    void init() {
      const FinalState fs(Cuts::abseta < 6.0);
      declare(fs, "FS");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.5), "Jets");

      VetoedFinalState fsv(fs);
      fsv.vetoNeutrinos();
      fsv.addVetoPair(PID::MUON);
      declare(fsv, "fsv");

      // For the MB ND selection
      const ChargedFinalState fschrgd(Cuts::abseta < 6.0);
      declare(fschrgd, "fschrgd");
      VetoedFinalState fschrgdv(fschrgd);
      fschrgdv.vetoNeutrinos();
      declare(fschrgdv, "fschrgdv");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        size_t offset = (en == "900"s)? 0 : 2;
        book(_h[en+"mb"],    1+offset, 1, 1); // energy flow in MB
        book(_h[en+"dijet"], 2+offset, 1, 1); // energy flow in dijet events
        book(_c[en+"mb"],    "/tmp/weightMB"+en);
        book(_c[en+"dijet"], "/tmp/weightDijet"+en);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }

    }


    void analyze(const Event& event) {
      // Skip if the event is empty
      const FinalState& fsv = apply<FinalState>(event, "fsv");
      if (fsv.empty()) vetoEvent;

      // Veto diffractive topologies according to defined hadron level
      double count_chrg_forward = 0;
      double count_chrg_backward = 0;
      const FinalState& fschrgdv = apply<FinalState>(event, "fschrgdv");
      for (const Particle& p : fschrgdv.particles()) {
        if (3.9 < p.eta() && p.eta() < 4.4)     ++count_chrg_forward;
        if (-4.4 < p.eta() && p.eta() < -3.9)  ++count_chrg_backward;
      }
      if (count_chrg_forward == 0 || count_chrg_backward == 0) vetoEvent;
      /// @todo "Diffractive" veto should really also veto dijet events?


      // MINIMUM BIAS EVENTS
      _c[_sqs+"mb"]->fill();
      for (const Particle& p: fsv.particles()) {
        _h[_sqs+"mb"]->fill(p.abseta(), p.E()/GeV);
      }


      // DIJET EVENTS
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const double PTCUT = (_sqs=="900"s)? 8*GeV : 20*GeV;
      const Jets jets = jetpro.jetsByPt(Cuts::pT > PTCUT);
      if (jets.size() >= 2) {
        // eta cut for the central jets
        if (jets[0].abseta() < 2.5 && jets[1].abseta() < 2.5) {
          // Back to back condition of the jets
          const double diffphi = deltaPhi(jets[1].phi(), jets[0].phi());
          if (diffphi-PI < 1.0) {
            _c[_sqs+"dijet"]->fill();
            for (const Particle& p: fsv.particles()) {
              _h[_sqs+"dijet"]->fill(p.abseta(), p.E()/GeV);
            }
          }
        }
      }

    }


    void finalize() {
      for (auto& item : _h) {
        scale(item.second, 0.5/_c[item.first]->sumW());
      }
    }


  private:

    /// @{
    map<string,Histo1DPtr> _h;
    map<string,CounterPtr> _c;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_I930319, CMS_2011_S9215166);

}
