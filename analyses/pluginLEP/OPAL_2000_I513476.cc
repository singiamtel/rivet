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
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(FastJets(fs, JetAlg::DURHAM, 0.7), "DurhamJets");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      if (isCompatibleWithSqrtS(172*GeV)) {
        ih = 0;
      } else if (isCompatibleWithSqrtS(183*GeV)) {
        ih = 1;
      } else if (isCompatibleWithSqrtS(189*GeV)) {
        ih = 2;
      }
      else {
        MSG_WARNING("Invalid CMS energy for OPAL_2000_I513476, assume 172 GeV.");
      }
      for (size_t i=0; i < _energies.size(); ++i) {
        book(_h[i]["thrust"]    ,  1,1,i+1);
        book(_h[i]["major"]     ,  2,1,i+1);
        book(_h[i]["minor"]     ,  3,1,i+1);
        book(_h[i]["aplanarity"],  4,1,i+1);
        book(_h[i]["oblateness"],  5,1,i+1);
        book(_h[i]["C"]         ,  6,1,i+1);
        book(_h[i]["rhoH"]      ,  7,1,i+1);
        book(_h[i]["sphericity"],  8,1,i+1);
        book(_h[i]["totalB"]    ,  9,1,i+1);
        book(_h[i]["wideB"]     , 10,1,i+1);
        book(_h[i]["y23"]       , 11,1,i+1);
        book(_h_mult[i]         , 13,1,i+1);
        book(_h[i]["pTin"]      , 15,1,i+1);
        book(_h[i]["pTout"]     , 16,1,i+1);
        book(_h[i]["y"]         , 17,1,i+1);
        book(_h[i]["x"]         , 18,1,i+1);
        book(_h[i]["xi"]        , 19,1,i+1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;

      // Thrust related
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _h[ih]["thrust"]    ->fill(thrust.thrust()     );
      _h[ih]["major"]     ->fill(thrust.thrustMajor());
      _h[ih]["minor"]     ->fill(thrust.thrustMinor());
      _h[ih]["oblateness"]->fill(thrust.oblateness() );

      // Sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h[ih]["sphericity"]->fill(sphericity.sphericity());
      _h[ih]["aplanarity"]->fill(sphericity.aplanarity());

      // C parameter
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h[ih]["C"]->fill(parisi.C());

      // Hemispheres
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");

      _h[ih]["rhoH"]  ->fill(hemi.scaledMhigh());
      _h[ih]["wideB"] ->fill(hemi.Bmax());
      _h[ih]["totalB"]->fill(hemi.Bsum());

      // Jets
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
      _h[ih]["y23"]->fill(y23);

      // charged particles
      _h_mult[ih]->fill(cfs.particles().size());
      for (const Particle& p : cfs.particles()) {
        const Vector3 mom3  = p.p3();
        const double energy = p.E();
        const double pTinT  = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
        _h[ih]["pTin"] ->fill(fabs(pTinT/GeV) );
        _h[ih]["pTout"]->fill(fabs(pToutT/GeV));
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
        _h[ih]["y"]->fill(fabs(rapidityT));
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        const double logInvScaledMom = -std::log(scaledMom);
        _h[ih]["xi"]->fill(logInvScaledMom);
        _h[ih]["x"] ->fill(scaledMom      );
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
      for (size_t i=0; i<_energies.size(); ++i) {
        const double sumw = _h_mult[i]->sumW();
        if (sumw == 0)  continue;

        scale(_h[i], 1./sumw);
        scale(_h_mult[i], 100./sumw);

        const double nch     = _h_mult[i]->xMean();
        const double nch_err = _h_mult[i]->xStdErr();
        m_ch->bin(i+1).set(nch, nch_err);
        double pTin     = _h[i]["pTin"]->xMean();
        double pTin_err = _h[i]["pTin"]->xStdErr();
        m_pTin->bin(i+1).set(pTin, pTin_err);
        double pTout     = _h[i]["pTout"]->xMean();
        double pTout_err = _h[i]["pTout"]->xStdErr();
        m_pTout->bin(i+1).set(pTout, pTout_err);
        double y     = _h[i]["y"]->xMean();
        double y_err = _h[i]["y"]->xStdErr();
        m_y->bin(i+1).set(y, y_err);
        double x     = _h[i]["x"]->xMean();
        double x_err = _h[i]["x"]->xStdErr();
        m_x->bin(i+1).set(x, x_err);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h[3];
    BinnedHistoPtr<int> _h_mult[3];
    size_t ih = 0;
    const vector<double> _energies{ 172., 183., 189. };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(OPAL_2000_I513476);


}
