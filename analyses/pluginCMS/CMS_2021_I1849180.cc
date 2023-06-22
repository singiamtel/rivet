// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Math/Vector3.hh"

namespace Rivet {

  /// @brief Drell-Yan dimuon production in proton-lead collisions at  8.16 TeV
  class CMS_2021_I1849180 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1849180);
    float y_shift = 0;

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs;

      const ParticlePair& beam = beams();
      pcom = beam.first.momentum() + beam.second.momentum();
      if (beam.first.mom().E() == beam.second.mom().E()) {
        y_shift = -0.465;
      }

      ZFinder zmumuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 10 * GeV, PID::MUON, 15.0 * GeV, 600.0 * GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zmumuFinder, "ZmumuFinder");

      ZFinder TotzmumuFinder(fs, Cuts::pT > 0.0 * GeV, PID::MUON, 0.0 * GeV, 1000.0 * GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(TotzmumuFinder, "TotzmumuFinder");

      declare(FinalState(), "FS");

      book(_h["111"], 1, 1, 1);
      book(_h["211"], 2, 1, 1);
      book(_h["311"], 3, 1, 1);
      book(_h["411"], 4, 1, 1);
      book(_h["511"], 5, 1, 1);
      book(_h["611"], 6, 1, 1);
      book(_h["711"], 7, 1, 1);
      book(_h["911"], 9, 1, 1);
      book(_h["1011"], 10, 1, 1);
      book(_h["1111"], 11, 1, 1);
      book(_h["1211"], 12, 1, 1);
      book(_h["1311"], 13, 1, 1);
      book(_h["1411"], 14, 1, 1);
      book(_h["1511"], 15, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ZFinder& ZmumuFinder = apply<ZFinder>(event, "ZmumuFinder");
      const ZFinder& TotzmumuFinder = apply<ZFinder>(event, "TotzmumuFinder");
      const Particles& zmumus = ZmumuFinder.bosons();
      const Particles& totzmumus = TotzmumuFinder.bosons();

      const Vector3 betacom = pcom.betaVec();
      const LorentzTransform comboost = LorentzTransform::mkFrameTransformFromBeta(betacom);

      if (zmumus.size() == 1) {
        for (const Particle& zmumu : zmumus) {
          Particles leptons = sortBy(ZmumuFinder.constituents(), cmpMomByPt);
          if (leptons[0].pt() < 15 or leptons[1].pt() < 10) {
            continue;
          }
          const FourMomentum lminus = leptons[0].charge() < 0 ? leptons[0].momentum() : leptons[1].momentum();
          const FourMomentum lplus = leptons[0].charge() < 0 ? leptons[1].momentum() : leptons[0].momentum();

          const double phi_acop = M_PI - deltaPhi(lminus, lplus);
          const double costhetastar = tanh((lminus.eta() - lplus.eta()) / 2.0);
          const double sin2thetastar = (costhetastar <= 1) ? 1.0 - sqr(costhetastar) : 0;
          const double phistar = tan(phi_acop / 2.0) * sqrt(sin2thetastar);

          const FourMomentum Zcm = comboost.transform(zmumu.momentum());

          const double Zmass = zmumu.mass() / GeV;
          const double Zpt = zmumu.momentum().pT() / GeV;

          const double Zy = Zcm.rapidity() + y_shift;

          if (Zcm.rapidity() < 1.93 && Zcm.rapidity() > -2.87) {
            _h["111"]->fill(Zmass);
          }
          if (Zmass < 60 && Zmass > 15) {
            if (Zcm.rapidity() < 1.93 && Zcm.rapidity() > -2.87) {
              _h["411"]->fill(Zpt);
              _h["611"]->fill(phistar);
            }
            _h["211"]->fill(Zy);
          }

          if (Zmass < 120 && Zmass > 60) {
            if (Zcm.rapidity() < 1.93 && Zcm.rapidity() > -2.87) {
              _h["511"]->fill(Zpt);
              _h["711"]->fill(phistar);
            }
            _h["311"]->fill(Zy);
          }
        }
      }

      if (totzmumus.size() == 1) {
        for (const Particle& totzmumu : totzmumus) {
          Particles leptons = sortBy(TotzmumuFinder.constituents(), cmpMomByPt);

          const FourMomentum lminus = leptons[0].charge() < 0 ? leptons[0].momentum() : leptons[1].momentum();
          const FourMomentum lplus = leptons[0].charge() < 0 ? leptons[1].momentum() : leptons[0].momentum();

          const double phi_acop = M_PI - deltaPhi(lminus, lplus);
          const double costhetastar = tanh((lminus.eta() - lplus.eta()) / 2.0);
          const double sin2thetastar = (costhetastar <= 1) ? 1.0 - sqr(costhetastar) : 0;
          const double totphistar = tan(phi_acop / 2.0) * sqrt(sin2thetastar);

          const FourMomentum totZcm = comboost.transform(totzmumu.momentum());

          const double totZmass = totzmumu.mass() / GeV;
          const double totZpt = totzmumu.momentum().pT() / GeV;
          const double totZy = totZcm.rapidity() + y_shift;

          if (totZcm.rapidity() < 1.93 && totZcm.rapidity() > -2.87) {
            _h["911"]->fill(totZmass);
          }

          if (totZmass < 60 && totZmass > 15) {
            if (totZcm.rapidity() < 1.93 && totZcm.rapidity() > -2.87) {
              _h["1211"]->fill(totZpt);
              _h["1411"]->fill(totphistar);
            }
            _h["1011"]->fill(totZy);
          }
          if (totZmass < 120 && totZmass > 60) {
            if (totZcm.rapidity() < 1.93 && totZcm.rapidity() > -2.87) {
              _h["1311"]->fill(totZpt);
              _h["1511"]->fill(totphistar);
            }
            _h["1111"]->fill(totZy);
          }
        }
      }
    }
    /// Normalise histograms etc., after the run
    void finalize() {
      double norm = crossSection() * 208 / nanobarn / sumW();

      scale(_h, norm);
    }

    FourMomentum pcom;
    map<string, Histo1DPtr> _h;
  };

  RIVET_DECLARE_PLUGIN(CMS_2021_I1849180);

}  // namespace Rivet
