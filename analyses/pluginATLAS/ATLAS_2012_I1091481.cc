// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Azimuthal ordering of charged hadrons
  class ATLAS_2012_I1091481 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2012_I1091481);

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs100(Cuts::abseta < 2.5 && Cuts::pT > 0.1*GeV);
      declare(cfs100,"CFS100");
      ChargedFinalState cfs500(Cuts::abseta < 2.5 && Cuts::pT > 0.5*GeV);
      declare(cfs500,"CFS500");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        size_t ih = bool(en == "900") + 1;

        book(_h[en+"E_10_100"], ih, 1, 1);
        book(_h[en+"E_1_100"],  ih, 1, 2);
        book(_h[en+"E_10_500"], ih, 1, 3);

        book(_h[en+"eta_10_100"], ih, 2, 1);
        book(_h[en+"eta_1_100"],  ih, 2, 2);
        book(_h[en+"eta_10_500"], ih, 2, 3);

        book(_c[en+"inclusive"], "_norm_inclusive_"+en);
        book(_c[en+"lowPt"],     "_norm_lowPt_"+en);
        book(_c[en+"pt500"],     "_norm_pt500_"+en);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    // Recalculate particle energy assuming pion mass
    double getPionEnergy(const Particle& p) {
      double m_pi = 0.1396*GeV;
      double p2 = p.p3().mod2()/(GeV*GeV);
      return sqrt(sqr(m_pi) + p2);
    }


    // S_eta core for one event
    //
    //  -1 + 1/Nch * |sum_j^Nch exp[i*(xi eta_j - Phi_j)]|^2
    //
    double getSeta(const Particles& part, double xi) {
      std::complex<double> c_eta (0.0, 0.0);
      for (const Particle& p : part) {
        double eta = p.eta();
        double phi = p.phi();
        double arg = xi*eta-phi;
         std::complex<double> temp(cos(arg), sin(arg));
         c_eta += temp;
      }
      return std::norm(c_eta)/part.size() - 1.0;
    }


    // S_E core for one event
    //
    //  -1 + 1/Nch * |sum_j^Nch exp[i*(omega X_j - Phi_j)]|^2
    //
    double getSE(const Particles& part, double omega) {
      double Xj = 0.0;
      std::complex<double> c_E (0.0, 0.0);
      for (unsigned int i=0; i < part.size(); ++i) {
        Xj += 0.5*getPionEnergy(part[i]);
        double phi = part[i].phi();
        double arg = omega*Xj - phi;
        std::complex<double> temp(cos(arg), sin(arg));
        c_E += temp;
        Xj += 0.5*getPionEnergy(part[i]);
      }
      return std::norm(c_E)/part.size() - 1.0;
    }


    // Convenient fill function
    void fillS(Histo1DPtr h, const Particles& part, bool SE=true) {
      // Loop over bins, take bin centers as parameter values
      // @todo use barchart()
      for (auto& b : h->bins()) {
        double x = b.xMid();
        double width = b.xWidth();
        double y;
        if (SE)  y = getSE(part,   x);
        else     y = getSeta(part, x);
        h->fill(x, y * width);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Charged fs
      const ChargedFinalState& cfs100  = apply<ChargedFinalState>(event, "CFS100");
      const Particles          part100 = cfs100.particles(cmpMomByEta);
      const ChargedFinalState& cfs500  = apply<ChargedFinalState>(event, "CFS500");
      const Particles&         part500 = cfs500.particles(cmpMomByEta);

      // Veto event if the most inclusive phase space has less than 10 particles and the max pT is > 10 GeV
      if (part100.size() < 11) vetoEvent;
      double ptmax = cfs100.particlesByPt()[0].pT()/GeV;
      if (ptmax > 10.0) vetoEvent;

      // Fill the pt>100, pTmax<10 GeV histos
      fillS(_h[_sqs+"E_10_100"],   part100, true);
      fillS(_h[_sqs+"eta_10_100"], part100, false);
      _c[_sqs+"inclusive"]->fill();

      // Fill the pt>100, pTmax<1 GeV histos
      if (ptmax < 1.0) {
        fillS(_h[_sqs+"E_1_100"],   part100, true);
        fillS(_h[_sqs+"eta_1_100"], part100, false);
        _c[_sqs+"lowPt"]->fill();
      }

      // Fill the pt>500, pTmax<10 GeV histos
      if (part500.size() > 10) {
        fillS(_h[_sqs+"E_10_500"],   part500, true );
        fillS(_h[_sqs+"eta_10_500"], part500, false);
        _c[_sqs+"pt500"]->fill();
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // The scaling takes the multiple fills per event into account
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));

        scale(_h[en+"E_10_100"], 1.0/ *_c[en+"inclusive"]);
        scale(_h[en+ "E_1_100"], 1.0/ *_c[en+"lowPt"]);
        scale(_h[en+"E_10_500"], 1.0/ *_c[en+"pt500"]);

        scale(_h[en+"eta_10_100"], 1.0/ *_c[en+"inclusive"]);
        scale(_h[en+ "eta_1_100"], 1.0/ *_c[en+"lowPt"]);
        scale(_h[en+"eta_10_500"], 1.0/ *_c[en+"pt500"]);
      }
    }


  private:

    map<string,Histo1DPtr> _h;
    map<string,CounterPtr> _c;

    string _sqs = "";
  };


  RIVET_DECLARE_PLUGIN(ATLAS_2012_I1091481);

}
