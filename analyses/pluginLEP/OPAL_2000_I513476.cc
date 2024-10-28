// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// @brief  event shapes at 172, 183, 189
  class OPAL_2000_I513476 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2000_I513476);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Projections
      const FinalState fs;
      declare(Beam(), "Beams");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(FastJets(fs, JetAlg::DURHAM, 0.7), "DurhamJets");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      size_t ih = 1;
      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal/MeV));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        book(_h[en+"thrust"],      1,1,ih);
        book(_h[en+"major"],       2,1,ih);
        book(_h[en+"minor"],       3,1,ih);
        book(_h[en+"aplanarity"],  4,1,ih);
        book(_h[en+"oblateness"],  5,1,ih);
        book(_h[en+"C"],           6,1,ih);
        book(_h[en+"rhoH"],        7,1,ih);
        book(_h[en+"sphericity"],  8,1,ih);
        book(_h[en+"totalB"],      9,1,ih);
        book(_h[en+"wideB"],      10,1,ih);
        book(_h[en+"y23"],        11,1,ih);
        book(_mult[en],           13,1,ih);
        book(_h[en+"pTin"],       15,1,ih);
        book(_h[en+"pTout"],      16,1,ih);
        book(_h[en+"y"],          17,1,ih);
        book(_h[en+"x"],          18,1,ih);
        book(_h[en+"xi"],         19,1,ih);
        ++ih;
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = 0.5*(beams.first.p3().mod() + beams.second.p3().mod());

      // Thrust related
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _h[_sqs+"thrust"]    ->fill(thrust.thrust()     );
      _h[_sqs+"major"]     ->fill(thrust.thrustMajor());
      _h[_sqs+"minor"]     ->fill(thrust.thrustMinor());
      _h[_sqs+"oblateness"]->fill(thrust.oblateness() );

      // Sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h[_sqs+"sphericity"]->fill(sphericity.sphericity());
      _h[_sqs+"aplanarity"]->fill(sphericity.aplanarity());

      // C parameter
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h[_sqs+"C"]->fill(parisi.C());

      // Hemispheres
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");

      _h[_sqs+"rhoH"]  ->fill(hemi.scaledMhigh());
      _h[_sqs+"wideB"] ->fill(hemi.Bmax());
      _h[_sqs+"totalB"]->fill(hemi.Bsum());

      // Jets
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
      _h[_sqs+"y23"]->fill(y23);

      // charged particles
      _mult[_sqs]->fill(cfs.particles().size());
      for (const Particle& p : cfs.particles()) {
        const Vector3 mom3  = p.p3();
        const double energy = p.E();
        const double pTinT  = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
        _h[_sqs+"pTin"] ->fill(fabs(pTinT/GeV) );
        _h[_sqs+"pTout"]->fill(fabs(pToutT/GeV));
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
        _h[_sqs+"y"]->fill(fabs(rapidityT));
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        const double logInvScaledMom = -std::log(scaledMom);
        _h[_sqs+"xi"]->fill(logInvScaledMom);
        _h[_sqs+"x"] ->fill(scaledMom);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // mean multiplicity
      BinnedEstimatePtr<int> m_ch;
      book(m_ch,14,1,1);
      // mean ptIn
      BinnedEstimatePtr<int> m_pTin;
      book(m_pTin,20,1,1);
      // mean ptOut
      BinnedEstimatePtr<int> m_pTout;
      book(m_pTout,20,1,2);
      // mean y
      BinnedEstimatePtr<int> m_y;
      book(m_y,20,1,3);
      // mean x
      BinnedEstimatePtr<int> m_x;
      book(m_x,20,1,4);

      // scale histos + fill averages
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal/MeV));

        const double sumw = _mult[en]->sumW();
        if (sumw == 0)  continue;
        scale(_mult[en], 100./sumw);

        for (auto& item : _h) {
          if (item.first.substr(0,6) != en)  continue;
          scale(item.second, 1./sumw);
        }

        for (size_t ih=1; ih <= m_ch->numBins(); ++ih) {
          if (m_ch->bin(ih).xEdge() != int(eVal))  continue;

          const double nch     = _mult[en]->xMean();
          const double nch_err = _mult[en]->xStdErr();
          m_ch->bin(ih).set(nch, nch_err);

          double pTin     = _h[en+"pTin"]->xMean();
          double pTin_err = _h[en+"pTin"]->xStdErr();
          m_pTin->bin(ih).set(pTin, pTin_err);

          double pTout     = _h[en+"pTout"]->xMean();
          double pTout_err = _h[en+"pTout"]->xStdErr();
          m_pTout->bin(ih).set(pTout, pTout_err);

          double y     = _h[en+"y"]->xMean();
          double y_err = _h[en+"y"]->xStdErr();
          m_y->bin(ih).set(y, y_err);

          double x     = _h[en+"x"]->xMean();
          double x_err = _h[en+"x"]->xStdErr();
          m_x->bin(ih).set(x, x_err);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h;

    map<string,BinnedHistoPtr<int>> _mult;

    string _sqs = "";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(OPAL_2000_I513476);


}
