// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"

namespace Rivet {


  /// @brief DELPHI event shapes below the Z pole
  class DELPHI_2003_I620250 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_2003_I620250);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections.
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");

      // Histogram booking offset numbers.
      unsigned int offset = 0;
      int offset2 = -1;
      isDisc = false;

      if      (isCompatibleWithSqrtS(45*GeV)) offset = 1;
      else if (isCompatibleWithSqrtS(66*GeV)) offset = 2;
      else if (isCompatibleWithSqrtS(76*GeV)) offset = 3;
      else if (isCompatibleWithSqrtS(183*GeV)) {
        offset2= 0;
        offset = 1;
      }
      else if (isCompatibleWithSqrtS(189*GeV)) {
        offset2= 0;
        offset = 2;
      }
      else if (isCompatibleWithSqrtS(192*GeV)) {
        offset2= 0;
        offset = 3;
      }
      else if (isCompatibleWithSqrtS(196*GeV)) {
        offset2= 0;
        offset = 4;
      }
      else if (isCompatibleWithSqrtS(200*GeV)) {
        offset2= 1;
        offset = 1;
      }
      else if (isCompatibleWithSqrtS(202*GeV)) {
        offset2= 1;
        offset = 2;
      }
      else if (isCompatibleWithSqrtS(205*GeV)) {
        offset2= 1;
        offset = 3;
      }
      else if (isCompatibleWithSqrtS(207*GeV)) {
        offset2= 1;
        offset = 4;
      }
      else    MSG_ERROR("Beam energy not supported!");
      // Book the histograms
      if (offset2 < 0) {
        book(_h["thrust"], 1, 1, offset);
        book(_h["major"], 2, 1, offset);
        book(_h["minor"], 3, 1, offset);
        book(_h["sphericity"], 4, 1, offset);
        book(_h["planarity"], 5, 1, offset);
        book(_h["oblateness"], 6, 1, offset);
        book(_h["heavy_jet_mass"], 7, 1, offset);
        book(_h["light_jet_mass"], 9, 1, offset);
        book(_h["diff_jet_mass"], 10, 1, offset);
        book(_h["total_jet_mass"], 11, 1, offset);
        book(_h["heavy_jet_mass_E"],  8, 1, offset);
        book(_h["total_jet_mass_E"], 12, 1, offset);
        book(_h["wide_broading"], 13, 1, offset);
        book(_h["narrow_broading"], 14, 1, offset);
        book(_h["total_broading"], 15, 1, offset);
        book(_h["diff_broading"], 16, 1, offset);
        book(_h["CParam"], 17, 1, offset);
      }
      else {
        isDisc = true;
        book(_d["rap"], 30+offset2, 1, offset);
        book(_d["xi"], 32+offset2, 1, offset);
        book(_d["pTIn"], 34+offset2, 1, offset);
        book(_d["pTOut"], 36+offset2, 1, offset);
        book(_d["thrust"], 38+offset2, 1, offset);
        book(_d["major"], 40+offset2, 1, offset);
        book(_d["minor"], 42+offset2, 1, offset);
        book(_d["oblateness"], 44+offset2, 1, offset);
        book(_d["wide_broading"], 46+offset2, 1, offset);
        book(_d["total_broading"], 48+offset2, 1, offset);
        book(_d["diff_broading"], 50+offset2, 1, offset);
        book(_d["CParam"], 52+offset2, 1, offset);
        book(_d["DParam"], 54+offset2, 1, offset);
        book(_d["heavy_jet_mass"], 56+offset2, 1, offset);
        book(_d["heavy_jet_mass_P"], 58+offset2, 1, offset);
        book(_d["heavy_jet_mass_E"], 60+offset2, 1, offset);
        book(_d["light_jet_mass"], 62+offset2, 1, offset);
        book(_d["diff_jet_mass"], 64+offset2, 1, offset);
        book(_d["sphericity"], 66+offset2, 1, offset);
        book(_d["planarity"], 68+offset2, 1, offset);
        book(_d["aplanarity"], 70+offset2, 1, offset);

        _axis["rap"] = YODA::Axis<double>({0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.5});
        _axis["xi"]  = YODA::Axis<double>({0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0,
                                           4.4, 4.8, 5.2, 5.6, 6.0, 6.4});
        _axis["pTIn"]  = YODA::Axis<double>({0.0, 0.1, 0.4, 0.65, 0.9, 1.1, 1.4, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0});
        _axis["pTOut"] = YODA::Axis<double>({0.0, 0.2, 0.4, 0.6, 0.85, 1.2, 1.6, 2.0, 3.0});
        _axis["thrust"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                              0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.24, 0.28, 0.32, 0.36});
        _axis["major"] = YODA::Axis<double>({0.0, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16,
                                             0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6});
        _axis["minor"] = YODA::Axis<double>({0.0, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.28, 0.32});
        _axis["oblateness"] = YODA::Axis<double>({0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16,
                                                  0.18, 0.20, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44});
        _axis["wide_broading"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                                                     0.08, 0.1, 0.12, 0.14, 0.17, 0.20, 0.24, 0.28});
        _axis["total_broading"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
                                                      0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.24, 0.27, 0.3, 0.33, 0.36});
        _axis["diff_broading"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                                                     0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.24, 0.28});
        _axis["CParam"] = YODA::Axis<double>({0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44,
                                              0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88});
        _axis["DParam"] = YODA::Axis<double>({0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44,
                                              0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92});
        _axis["heavy_jet_mass"] = YODA::Axis<double>({0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16,
                                                      0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.54});
        _axis["heavy_jet_mass_P"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08,
                                                        0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.28, 0.32});
        _axis["heavy_jet_mass_E"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05});
        _axis["light_jet_mass"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05});
        _axis["diff_jet_mass"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3});
        _axis["sphericity"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1,
                                                  0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6});
        _axis["planarity"] = YODA::Axis<double>({0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08,
                                                 0.1, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4});
        _axis["aplanarity"] = YODA::Axis<double>({0.0, 0.004, 0.01, 0.016, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.12, 0.16});
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        for (const auto& item : _axis) {
          _edges[item.first] = _d[item.first]->xEdges();
        }
      }

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      smartfill("thrust", 1.-thrust.thrust());
      smartfill("major", thrust.thrustMajor());
      smartfill("minor", thrust.thrustMinor());
      smartfill("oblateness", thrust.oblateness() );
      smartfill("sphericity", sphericity.sphericity());
      smartfill("planarity", sphericity.planarity() );
      smartfill("aplanarity", sphericity.aplanarity());
      smartfill("heavy_jet_mass", hemi.scaledM2high());
      smartfill("light_jet_mass", hemi.scaledM2low() );
      smartfill("diff_jet_mass", hemi.scaledM2diff());
      smartfill("wide_broading", hemi.Bmax());
      smartfill("narrow_broading", hemi.Bmin());
      smartfill("total_broading", hemi.Bsum());
      smartfill("diff_broading", hemi.Bdiff());
      smartfill("total_jet_mass", hemi.scaledM2low()+hemi.scaledM2high());
      // E and p scheme jet masses
      Vector3 axis = thrust.thrustAxis();
      FourMomentum p4WithE, p4AgainstE;
      FourMomentum p4WithP, p4AgainstP;
      double Evis(0);
      for (const Particle& p : apply<FinalState>(event, "FS").particles()) {
        Vector3 p3 = p.momentum().vector3().unitVec();
        const double   E = p.momentum().E();
        Evis += E;
        p3 = E*p3;
        const double p3Para = dot(p3, axis);
        FourMomentum p4E(E,p3.x(),p3.y(),p3.z());
        FourMomentum p4P(p.p3().mod(),p.p3().x(),p.p3().y(),p.p3().z());
        if (p3Para > 0)      {
          p4WithE    += p4E;
          p4WithP    += p4P;
        }
        else if (p3Para < 0) {
          p4AgainstE += p4E;
          p4AgainstP += p4P;
        }
        else {
          MSG_WARNING("Particle split between hemispheres");
          p4WithE    += 0.5 * p4E;
          p4AgainstE += 0.5 * p4E;
          p4WithP    += 0.5 * p4P;
          p4AgainstP += 0.5 * p4P;
        }
      }
      // E scheme
      const double mass2With_E    = p4WithE.mass2()/sqr(Evis);
      const double mass2Against_E = p4AgainstE.mass2()/sqr(Evis);
      // fill the histograms
      smartfill("heavy_jet_mass_E", max(mass2With_E,mass2Against_E));
      smartfill("total_jet_mass_E", mass2With_E+mass2Against_E);
      // pscheme
      const double mass2With_P    = p4WithP.mass2()/sqr(Evis);
      const double mass2Against_P = p4AgainstP.mass2()/sqr(Evis);
      // fill the histograms
      smartfill("heavy_jet_mass_P", max(mass2With_P, mass2Against_P));

      MSG_DEBUG("Calculating Parisi params");
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      smartfill("CParam", parisi.C());
      smartfill("DParam", parisi.D());

      // single particle distributions
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (isDisc) {
        for (const Particle& p : fs.particles()) {
          if ( ! PID::isCharged(p.pid())) continue;
          // Get momentum and energy of each particle.
          const Vector3 mom3 = p.p3();
          const double energy = p.E();

          // Scaled momenta.
          const double mom = mom3.mod();
          const double scaledMom = mom/meanBeamMom;
          const double logInvScaledMom = -std::log(scaledMom);
          smartfill("xi", logInvScaledMom);

          // Get momenta components w.r.t. thrust and sphericity.
          const double momT = dot(thrust.thrustAxis(), mom3);
          const double pTinT = dot(mom3, thrust.thrustMajorAxis());
          const double pToutT = dot(mom3, thrust.thrustMinorAxis());
          smartfill("pTIn", fabs(pTinT/GeV));
          smartfill("pTOut", fabs(pToutT/GeV));

          // Calculate rapidities w.r.t. thrust and sphericity.
          const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
          smartfill("rap", fabs(rapidityT));
          MSG_TRACE(fabs(rapidityT) << " " << scaledMom/GeV);
        }
      }
    }

    void smartfill(const string& tag, const double value) {
      if (isDisc) {
        const size_t idx = _axis[tag].index(value);
        if (idx && idx <= _edges[tag].size()) {
          _d[tag]->fill(_edges[tag][idx-1]);
        }
        else {
          _d[tag]->fill(string("OTHER"));
        }
      }
      else {
        _h[tag]->fill(value);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
      for (auto& item : _d) {
        if (item.first == "rap" ||
            item.first == "xi"  ||
            item.first == "pTIn" ||
            item.first == "pTOut")  scale(item.second, 1./sumOfWeights());
        else   normalize(item.second);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, BinnedHistoPtr<string>> _d;
    map<string, YODA::Axis<double>> _axis;
    map<string, vector<string>> _edges;
    bool isDisc;

    /// @}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_2003_I620250);


}
