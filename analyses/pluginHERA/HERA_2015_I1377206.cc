// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Measurement of sigma_red (F2) of H1 and ZEUS at different beam energies
  class HERA_2015_I1377206 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HERA_2015_I1377206);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");

      Histo1DPtr dummy;
      string beamOpt = getOption<string>("BEAM","NONE");

      if (beamOpt == "NONE") {
        const ParticlePair& beam = beams();
        _positron = (beam.first.pid() == PID::POSITRON || beam.second.pid() == PID::POSITRON );
      } else {
        if (beamOpt == "EMINUS") _positron = false;
        else if (beamOpt == "EPLUS") _positron = true;
        else {
          MSG_ERROR("Beam option error. You have specified an unsupported beam.");
          return;
        }
      }

      // Book beams-dependent histograms
      const double eps = 0.01;
      if (isCompatibleWithSqrtS(318*GeV, eps) && _positron  ) {
        // NC e+ p at sqrts=318
        const vector<double> Q2edges = {
          0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.6, 0.7, 1.0, 1.1, 1.3, 1.7, 2.3,
          3.1, 3.8, 5.3, 8.0, 9.1, 11., 13., 17.4, 19.1, 25.8, 28., 30., 42.,
          49., 54., 65., 75., 108., 134., 180., 225., 280., 325., 355., 455.,460.,
          545., 560., 765., 770., 835., 900., 1120., 1295., 1300., 1755., 1800.,
          2270., 2500., 3685., 4000., 6520., 7000., 9275., 10000., 15000., 17000.,
          24770., 25000., 42000.
        };
        book(_h_sigred, Q2edges);
        _h_sigred->maskBins({9, 24, 27, 36, 38, 39, 41, 44, 46, 48, 50, 52, 54, 56, 58});
        size_t idx = 0;
        for (auto& b : _h_sigred->bins()) {
          book(b, 1, 1, ++idx);
        }
        // CC e+ p at sqrts=318
        book(_h_sigred_cc, {280., 325., 460., 545., 900., 1120., 1300., 1755., 1800., 2270.,
                            2500., 3685., 4000., 6520., 7000., 9275., 10000., 20000., 42000.});
        _h_sigred_cc->maskBins({2, 4, 6, 8, 10, 12, 14, 16});
        idx = 0;
        for (auto& b : _h_sigred_cc->bins()) {
          book(b, 6, 1, ++idx);
        }
      }
      else if (isCompatibleWithSqrtS(300*GeV, eps) && _positron) {
        // NC e+ p at sqrts=300
        const vector<double> Q2edges = {
          0.01, 0.05, 0.07, 0.09, 0.12, 0.18, 0.22, 0.32, 0.4, 0.45, 0.6, 0.7, 1.0, 1.1, 1.3,
          1.7, 2.3, 3.1, 3.8, 5.3, 8., 9.1, 11., 13., 17.4, 19.1, 25.8, 28., 30., 42., 49.,
          54., 65., 75., 108., 134., 180., 225., 280., 325., 355., 455., 460., 545., 560.,
          765., 770., 835., 900., 1120., 1295., 1300., 1755., 1800., 2270., 2500., 3685.,
          4000., 6520., 7000., 9275., 10000., 15000., 17000., 24770., 25000., 42000.
        };
        book(_h_sigred, Q2edges);
        _h_sigred->maskBins({13, 28, 31, 40, 42, 44, 46, 48, 51, 53, 55, 57, 59, 61, 63, 65});
        size_t idx = 0;
        for (auto& b : _h_sigred->bins()) {
          book(b, 2, 1, ++idx);
        }
      }
      else if (isCompatibleWithSqrtS(251*GeV, eps) && _positron  ) {
        // NC e+ p at sqrts=251
        const vector<double> Q2edges = {
          1., 1.7, 2.3, 3.1, 3.8, 5.3, 8., 9.1, 11., 13., 17.4, 22.1, 28.,
          30., 42., 49., 54., 65., 75., 108., 134., 180., 225., 280., 325.,
          355., 455., 460., 545., 560., 765., 770., 835.
        };
        book(_h_sigred, Q2edges);
        _h_sigred->maskBins({8, 13, 16, 18, 25, 29, 31});
        size_t idx = 0;
        for (auto& b : _h_sigred->bins()) {
          book(b, 3, 1, ++idx);
        }
      }

      // e+- p at sqrts=225
      else if (isCompatibleWithSqrtS(225*GeV, eps)) {
        if (_positron) {
          // NC e+ p at sqrts=225
          const vector<double> Q2edges = {
            1., 1.7, 2.3, 3.1, 3.8, 5.3, 8., 9.1, 11., 13., 17.4, 22.1, 28.,
            30., 42., 49., 54., 65., 75., 108., 134., 180., 225., 280., 325.,
            355., 455., 460., 545., 560., 765., 770., 835.
          };
          book(_h_sigred, Q2edges);
          _h_sigred->maskBins({8, 13, 16, 18, 25, 29, 31});
          size_t idx = 0;
          for (auto& b : _h_sigred->bins()) {
            book(b, 4, 1, ++idx);
          }
        }
        else {
          // NC e- p at sqrts=225
          const vector<double> Q2edges = {
            54., 65., 75., 108., 134., 180., 225., 280., 325., 355., 455.,
            460., 545., 560., 765., 770., 835., 900., 1120., 1295., 1300.,
            1755., 1800., 2270., 2500., 3685., 4000., 6520., 7000., 9275.,
            10000., 15000., 17000., 24770., 25000., 42000., 70000.
          };
          book(_h_sigred, Q2edges);
          _h_sigred->maskBins({2, 9, 11, 13, 15, 17, 20, 22, 24, 26, 28, 30, 32, 34});
          size_t idx = 0;
          for (auto& b : _h_sigred->bins()) {
            book(b, 5, 1, ++idx);
          }
          // CC e- p at sqrts=225
          book(_h_sigred_cc, {280., 325., 460., 545., 900., 1120., 1300., 1755., 1800., 2270.,
                              2500., 3685., 4000., 6520., 7000., 9275., 10000., 20000., 42000.});
          _h_sigred_cc->maskBins({2, 4, 6, 8, 10, 12, 14, 16});
          idx = 0;
          for (auto& b : _h_sigred_cc->bins()) {
            book(b, 7, 1, ++idx);
          }
        }
      }
    }


    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;

      // Flux factor
      const double alpha = 7.29927e-3;
      // GF = 1.16638e-5 Fermi constant
      const double GF2 = 1.16638e-5*1.16638e-5;
      // MW = 80.385 W-boson mass
      const double MW2 = 80.385 * 80.385;

      if (PID::isNeutrino(dl.out().abspid()) ) {
        // fill histo for CC
        double F = 2.0*M_PI*x/GF2 * sqr((MW2 + Q2)/MW2);
        _h_sigred_cc->fill(Q2,x,F); // fill histogram x,Q2
      }
      else {
        // fill histo for NC
        double F = x*sqr(Q2)/(2.0*M_PI*sqr(alpha)*(1.0+sqr(1-y)));
        _h_sigred->fill(Q2,x,F); // fill histogram x,Q2
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double gev2nb = 0.389e6;
      const double scalefactor=crossSection()/nanobarn/sumOfWeights()/gev2nb ;
      // with _h_sigred.scale also q2 bin width is scaled
      scale(_h_sigred, scalefactor);
      scale(_h_sigred_cc, scalefactor);
      divByGroupWidth({_h_sigred, _h_sigred_cc});
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_sigred, _h_sigred_cc;
    Histo1DPtr _hist_Q2_10,_hist_Q2_100,_hist_Q2_1000,_hist_Q2_2000,_hist_Q2_3000;
    bool _positron;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(HERA_2015_I1377206);

}
