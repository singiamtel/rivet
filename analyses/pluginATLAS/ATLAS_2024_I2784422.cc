// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {

  /// @brief ATLAS Kaon and Lambda Underlying Event spectra at 13 TeV
  class ATLAS_2024_I2784422 : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2024_I2784422);

    void init() {

      // Jet finding over all ATLAS inner detector prompt-charged above pT threshold.
      declare(FastJets(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV), JetAlg::ANTIKT, 0.4), "jets");

      // We will only then use prompt charged at a smaller radius than the jet maximum: 2.5 - R
      // This is so as to not be biased by charged particles from any jets at larger radii.
      declare(ChargedFinalState(Cuts::abseta < 2.1 && Cuts::pT > 500*MeV), "promptCharged");

      // But we also have an event-requirement of at least one 1 GeV track somewhere in the inner detector
      declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 1*GeV), "eventSelection");

      // We take strange particles at a tighter eta still, this was to preserve a high purity from the detector reconstruction.
      declare(UnstableFinalState(Cuts::pid == PID::K0S && Cuts::abseta < 1.0), "kaon");

      declare(UnstableFinalState(Cuts::abspid == PID::LAMBDA && Cuts::abseta < 1.0), "lambda");

      // All 1D observables used in the creation of final ratios.
      size_t ih = 1;
      for (const string& obs : vector<string>{ "pTlead", "nTrans" }) {
        for (const string& num : vector<string>{ "kaon", "lambda" }) {
          for (const string& den : vector<string>{ "", "prompt", "kaon", "lambda" }) {
            for (const string& var : variables) {
              for (const string& reg : regions) {
                if (den == "")  dualbook(obs, num, var, reg, ih);
                else            dualbook(obs, num, den, var, reg, ih);
                ih += 2;
              }
              if (num == den)  break; // nTrans only
            }
            if (num == den)  break;
          }
        }
      }
    }

    void analyze(const Event& e) {

      // Require at least one charged particle to be above 1 GeV
      const Particles& eventSelection = apply<ChargedFinalState>(e, "eventSelection").particles();
      if (eventSelection.empty())  vetoEvent;

      // We want our leading jet to be fully contained within the ATLAS inner detector, hence cutting jets at eta=2.1 for R=0.4
      const Jets& jets = apply<FastJets>(e, "jets").jetsByPt(Cuts::absrap < 2.1 && Cuts::pT > 1*GeV);
      if (jets.empty())  vetoEvent;

      const double jetphi = jets[0].phi();
      const double jetpT  = jets[0].pT();

      // Count the multiplicity and scalar sum-pt of charged particles and strange particles
      map<string,YODA::Counter> counts;

      // Charged particles up to the max radius allowed for the leading jet
      const Particles& promptCharged = apply<ChargedFinalState>(e, "promptCharged").particles();
      for (const Particle& p : promptCharged) {
        if (p.origin().polarRadius() / mm > 1.5) continue;
        // Counting prompties
        const double dPhi = deltaPhi(p.phi(), jetphi);
        if (dPhi <= PI/3.0) {
          counts["prompt"s+"towards"s+"n"s].fill();
          counts["prompt"s+"towards"s+"sumpt"s].fill(p.pT()/GeV);
        } else if (dPhi > 2*PI/3.0) {
          counts["prompt"s+"away"s+"n"s].fill();
          counts["prompt"s+"away"s+"sumpt"s].fill(p.pT()/GeV);
        } else {
          counts["prompt"s+"transverse"s+"n"s].fill();
          counts["prompt"s+"transverse"s+"sumpt"s].fill(p.pT()/GeV);
        }
      }

      for (const string& type : vector<string>{"lambda", "kaon"}) {

        // Strange particles up to eta of 1.0
        const Particles& identifiedFS = apply<UnstableFinalState>(e, type).particles();

        const double RxyMin = (type == "kaon"s? 4*mm : 17*mm);
        const double pTMin  = (type == "kaon"s? 400*MeV : 750*MeV);

        // Loop over Kaon and Lambda particle vectors, apply detector cuts
        for (const Particle& p : identifiedFS) {

          const double particlePt = p.pT();
          if (particlePt < pTMin)  continue;

          const Particles children = p.children();
          if (children.size() != 2)  continue;

          if (children[0].abseta() > 2.5 || children[1].abseta() > 2.5)  continue;

          const double decayDistancePerp = children[0].origin().polarRadius() / mm;
          if (decayDistancePerp <= RxyMin || decayDistancePerp > 300*mm)  continue;

          if (type == "kaon"s && (children[0].abspid() != PID::PIPLUS || children[1].abspid() != PID::PIPLUS)) {
            continue;
          }
          bool lambdaRejected = false;
          if (children[0].abspid() != PID::PIPLUS && children[0].abspid() != PID::PROTON) lambdaRejected = true;
          if (children[1].abspid() != PID::PIPLUS && children[1].abspid() != PID::PROTON) lambdaRejected = true;
          if (type == "lambda"s && lambdaRejected)  continue;


          const double dPhi = deltaPhi(p.phi(), jetphi);
          if (dPhi < PI/3.0) {
            counts[type+"towards"+"n"].fill();
            counts[type+"towards"+"sumpt"].fill(p.pT()/GeV);
          } else if (dPhi > 2*PI/3.0) {
            counts[type+"away"+"n"].fill();
            counts[type+"away"+"sumpt"].fill(p.pT()/GeV);
          } else {
            counts[type+"transverse"+"n"].fill();
            counts[type+"transverse"+"sumpt"].fill(p.pT()/GeV);
          }
        }

      }

      const bool isRestrictedRange = (jetpT > 10*GeV && jetpT <= 40*GeV);
      const double nTransversePrompt = counts["prompt"s+"transverse"s+"n"s].val();

      // Fill all 1D accumulator histograms
      for (const string& type : vector<string>{ "kaon", "lambda", "prompt" }) {
        for (const string& var : variables) {
          for (const string& reg : regions) {
            _h["pTlead"+type+reg+var]->fill(jetpT/GeV, counts[type+reg+var].val());
            if (isRestrictedRange) {
              _h["nTrans"+type+reg+var]->fill(nTransversePrompt, counts[type+reg+var].val());
            }
          }
        }
      }

      // Fill the two "per event" normalisation histograms
      _h["pTlead"]->fill(jetpT/GeV);
      if (isRestrictedRange)  _h["nTrans"]->fill(nTransversePrompt);

    }

    void finalize() {

      for (const string& obs : vector<string>{ "pTlead", "nTrans" }) {
        for (const string& num : vector<string>{ "kaon", "lambda" }) {
          for (const string& reg : regions) {
            // Normalise same-species ratio distributions,
            // used to obtain ensemble-averaged mean-pT
            divide(_h[obs+num+reg+"sumpt"], _h[obs+num+reg+"n"], _e[obs+num+num+reg+"n"]);

            // Normalise inter-species ratio distributions
            for (const string& den : vector<string>{ "prompt", "kaon" }) {
              if (num == den)  break;
              for (const string& var : variables) {
                divide(_h[obs+num+reg+var], _h[obs+den+reg+var], _e[obs+num+den+reg+var]);
              }
            }
          } // regions
        }

        // Normalise per-event normalised distributions
        // Note: Due to scale(), this must run AFTER the inter-species ratios above,
        // and that kPrompt has a different normalisation factor to kKaon or kLambda.
        for (const string& num : vector<string>{ "kaon", "lambda" }) {
          for (const string& reg : regions) {
            const double deltaEta = 2.*1.0;
            const double deltaPhi = 2.*M_PI/3.;
            const double regionSize = deltaEta * deltaPhi;
            for (const string& var : variables) {
              scale(_h[obs+num+reg+var], 1.0/regionSize);
              divide(_h[obs+num+reg+var], _h[obs], _e[obs+num+reg+var]);
            }
          }
        }

      } // pTlead vs nTrans

    }

    void dualbook(const string& obs, const string& num, const string& var, const string& reg, size_t ih) {
      const string label = obs+num+reg+var;
      book(_e[label], ih, 1, 1);
      if (_h.find(obs) == _h.end()) {
        book(_h[obs], "_"+obs+"Spectrum", _e[label].binning().edges<0>());
      }
    }

    void dualbook(const string& obs, const string& num, const string& den, const string& var, const string& reg, size_t ih) {
      const string ratio = obs+num+den+reg+var;
      const string numer = obs+num+reg+var;
      const string denom = obs+den+reg+var;
      book(_e[ratio], ih, 1, 1);
      if (_h.find(numer) == _h.end())  book(_h[numer], "_"+obs+"_"+num+"_"+reg+"_"+var, _e[ratio].binning().edges<0>());
      if (_h.find(denom) == _h.end())  book(_h[denom], "_"+obs+"_"+den+"_"+reg+"_"+var, _e[ratio].binning().edges<0>());
    }

    private:

    const vector<string> regions = {"away", "towards", "transverse"};
    const vector<string> variables = {"n", "sumpt"};

    map<string,Histo1DPtr> _h;
    map<string,Estimate1DPtr> _e;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2024_I2784422);

}
