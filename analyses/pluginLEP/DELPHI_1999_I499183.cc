// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"

namespace Rivet {


  /// @brief event shapes at 133, 161 172 and 183 GeV
  class DELPHI_1999_I499183 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1999_I499183);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      size_t ie=0;
      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        size_t offset , offset2;
        if (ie==3)      { offset = 1; offset2 = 1; }
        else if (ie==2) { offset = 0; offset2 = 3; }
        else if (ie==1) { offset = 0; offset2 = 2; }
        else            { offset = 0; offset2 = 1; }

        book(_h[en+"thrust"]         , 13+offset, 1, offset2);
        book(_h[en+"major"]          , 15+offset, 1, offset2);
        book(_h[en+"minor"]          , 17+offset, 1, offset2);
        book(_h[en+"oblateness"]     , 19+offset, 1, offset2);
        book(_h[en+"sphericity"]     , 21+offset, 1, offset2);
        book(_h[en+"planarity"]      , 23+offset, 1, offset2);
        book(_h[en+"aplanarity"]     , 25+offset, 1, offset2);
        book(_h[en+"heavy_jet_mass"] , 27+offset, 1, offset2);
        book(_h[en+"light_jet_mass"] , 29+offset, 1, offset2);
        book(_h[en+"diff_jet_mass"]  , 31+offset, 1, offset2);
        book(_h[en+"wide_broading"]  , 33+offset, 1, offset2);
        book(_h[en+"narrow_broading"], 35+offset, 1, offset2);
        book(_h[en+"total_broading"] , 37+offset, 1, offset2);
        book(_h[en+"diff_broading"]  , 39+offset, 1, offset2);
        book(_h[en+"CParam"]         , 41+offset, 1, offset2);
        book(_h[en+"DParam"]         , 43+offset, 1, offset2);

        if (ie < 3) {
          book(_p["thrust"][ie], 1,1,1+ie);
          book(_p["major"][ie],  2,1,1+ie);
          book(_p["minor"][ie],  3,1,1+ie);
          book(_p["obl"][ie],    4,1,1+ie);
          book(_p["heavy"][ie],  5,1,1+ie);
          book(_p["light"][ie],  6,1,1+ie);
          book(_p["diff"][ie],   7,1,1+ie);
          book(_p["bmax"][ie],   8,1,1+ie);
          book(_p["bmin"][ie],   9,1,1+ie);
          book(_p["bsum"][ie],  10,1,1+ie);
          book(_p["bdiff"][ie], 11,1,1+ie);
          book(_p["C"][ie],     12,1,1+ie);
        }
        ++ie;
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }

    void fillMoment(array<BinnedProfilePtr<int>,3>& mom, double val) {
      double tmp=val;
      for (size_t ix=0; ix<3; ++ix) {
        mom[ix]->fill(std::stoi(_sqs), tmp);
        tmp *= val;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = 0.5*(beams.first.p3().mod() + beams.second.p3().mod());
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // thrust related observables
      _h[_sqs+"thrust"]    ->fill(1.-thrust.thrust()  );
      _h[_sqs+"major"]     ->fill(thrust.thrustMajor());
      _h[_sqs+"minor"]     ->fill(thrust.thrustMinor());
      _h[_sqs+"oblateness"]->fill(thrust.oblateness() );
      fillMoment(_p["thrust"],1.-thrust.thrust()  );
      fillMoment(_p["major"] ,thrust.thrustMajor());
      fillMoment(_p["minor"] ,thrust.thrustMinor());
      fillMoment(_p["obl"]   ,thrust.oblateness() );
      // sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h[_sqs+"sphericity"]->fill(sphericity.sphericity());
      _h[_sqs+"planarity"] ->fill(sphericity.planarity() );
      _h[_sqs+"aplanarity"]->fill(sphericity.aplanarity());
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      // standard jet masses
      _h[_sqs+"heavy_jet_mass"]->fill(hemi.scaledM2high());
      _h[_sqs+"light_jet_mass"]->fill(hemi.scaledM2low() );
      _h[_sqs+"diff_jet_mass"] ->fill(hemi.scaledM2diff());
      fillMoment(_p["heavy"],hemi.scaledM2high());
      fillMoment(_p["light"],hemi.scaledM2low() );
      fillMoment(_p["diff"] ,hemi.scaledM2diff());
      // jet broadening
      _h[_sqs+"wide_broading"]  ->fill(hemi.Bmax() );
      _h[_sqs+"narrow_broading"]->fill(hemi.Bmin() );
      _h[_sqs+"total_broading"] ->fill(hemi.Bsum() );
      _h[_sqs+"diff_broading"]  ->fill(hemi.Bdiff());
      fillMoment(_p["bmax"],hemi.Bmax() );
      fillMoment(_p["bmin"],hemi.Bmin() );
      fillMoment(_p["bsum"],hemi.Bsum() );
      fillMoment(_p["bdiff"],hemi.Bdiff());
      MSG_DEBUG("Calculating Parisi params");
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h[_sqs+"CParam"]->fill(parisi.C());
      _h[_sqs+"DParam"]->fill(parisi.D());
      fillMoment(_p["C"],parisi.C());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));

        for (auto& item : _h) {
          if (item.first.substr(0,3) != en)  continue;
          normalize(item.second);
        }

      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h;

    map<string,array<BinnedProfilePtr<int>,3>> _p;

    string _sqs = "";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_1999_I499183);


}
