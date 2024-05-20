// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Charged particle multiplicities in deep inelastic scattering at HERA (H1)
  class H1_1996_I422230 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1996_I422230);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const DISLepton disl;
      declare(disl, "Lepton");
      declare(DISKinematics(), "Kinematics");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs;
      declare(fs,"FS");
      const ChargedFinalState cfs(disl.remainingFinalState());
      declare(cfs,"CFS");

      //binned histograms to count for the multiplicity
      for (size_t ix = 0; ix < 4; ++ix) {
        book(_Nevt_after_cuts[ix], "TMP/Nevt_after_cuts"+ to_string(ix));
      }

      const vector<double> WEdges{80., 115, 150., 185., 220.};
      book(_g["mult1"], WEdges);
      book(_g["mult2"], WEdges);
      book(_g["mult3"], WEdges);
      book(_g["mult4"], WEdges);
      book(_g["mult_all"], WEdges);
      book(_g["mult10_all"], WEdges);
      book(_g["mult11_all"], WEdges);
      book(_g["mult12_all"], WEdges);
      for (size_t iW=0; iW < _g["mult1"]->numBins(); ++iW) {
        book(_g["mult1"]->bin(iW+1), iW+1, 1, 1);
        book(_g["mult2"]->bin(iW+1), iW+1, 1, 2);
        book(_g["mult3"]->bin(iW+1), iW+1, 1, 3);
        book(_g["mult4"]->bin(iW+1), iW+1, 1, 4);
        const auto& ref = refData<YODA::BinnedEstimate<int>>(4,1,4);
        book(_g["mult_all"]->bin(iW+1),   "TMP/dummy" + to_string(iW), ref);
        book(_g["mult10_all"]->bin(iW+1), "TMP/dummy1"+ to_string(iW), ref);
        book(_g["mult11_all"]->bin(iW+1), "TMP/dummy2"+ to_string(iW), ref);
        book(_g["mult12_all"]->bin(iW+1), "TMP/dummy3"+ to_string(iW), ref);
      }

      //histograms for the statistical moments and the mean
      book(_e["mean0"],5,1,1);
      book(_e["D2_0"],5,1,2);
      book(_e["D3_0"],5,1,3);
      book(_e["D4_0"],5,1,4);
      book(_e["C2_0"],5,1,5);
      book(_e["C3_0"],5,1,6);
      book(_e["C4_0"],5,1,7);
      book(_e["R2_0"],5,1,8);
      book(_e["R3_0"],5,1,9);

      book(_e["mean12"],6,1,1);
      book(_e["D2_12"],6,1,2);
      book(_e["D3_12"],6,1,3);
      book(_e["D4_12"],6,1,4);
      book(_e["C2_12"],6,1,5);
      book(_e["C3_12"],6,1,6);
      book(_e["C4_12"],6,1,7);
      book(_e["R2_12"],6,1,8);
      book(_e["R3_12"],6,1,9);
      book(_e["K3_12"],6,1,10);

      book(_e["mean13"],7,1,1);
      book(_e["D2_13"],7,1,2);
      book(_e["D3_13"],7,1,3);
      book(_e["D4_13"],7,1,4);
      book(_e["C2_13"],7,1,5);
      book(_e["C3_13"],7,1,6);
      book(_e["C4_13"],7,1,7);
      book(_e["R2_13"],7,1,8);
      book(_e["R3_13"],7,1,9);
      book(_e["K3_13"],7,1,10);

      book(_e["mean14"],8,1,1);
      book(_e["D2_14"],8,1,2);
      book(_e["D3_14"],8,1,3);
      book(_e["D4_14"],8,1,4);
      book(_e["C2_14"],8,1,5);
      book(_e["C3_14"],8,1,6);
      book(_e["C4_14"],8,1,7);
      book(_e["R2_14"],8,1,8);
      book(_e["R3_14"],8,1,9);

      book(_e["mean15"],9,1,1);
      book(_e["D2_15"],9,1,2);
      book(_e["D3_15"],9,1,3);
      book(_e["D4_15"],9,1,4);
      book(_e["C2_15"],9,1,5);
      book(_e["C3_15"],9,1,6);
      book(_e["C4_15"],9,1,7);
      book(_e["R2_15"],9,1,8);
      book(_e["R3_15"],9,1,9);

      book(_e["mean23"],10,1,1);
      book(_e["D2_23"],10,1,2);
      book(_e["D3_23"],10,1,3);
      book(_e["D4_23"],10,1,4);
      book(_e["C2_23"],10,1,5);
      book(_e["C3_23"],10,1,6);
      book(_e["C4_23"],10,1,7);
      book(_e["R2_23"],10,1,8);
      book(_e["R3_23"],10,1,9);
      book(_e["K3_23"],10,1,10);

      book(_e["mean34"],11,1,1);
      book(_e["D2_34"],11,1,2);
      book(_e["D3_34"],11,1,3);
      book(_e["D4_34"],11,1,4);
      book(_e["C2_34"],11,1,5);
      book(_e["C3_34"],11,1,6);
      book(_e["C4_34"],11,1,7);
      book(_e["R2_34"],11,1,8);
      book(_e["R3_34"],11,1,9);
      book(_e["K3_34"],11,1,10);

      book(_e["mean45"],12,1,1);
      book(_e["D2_45"],12,1,2);
      book(_e["D3_45"],12,1,3);
      book(_e["D4_45"],12,1,4);
      book(_e["C2_45"],12,1,5);
      book(_e["C3_45"],12,1,6);
      book(_e["C4_45"],12,1,7);
      book(_e["R2_45"],12,1,8);
      book(_e["R3_45"],12,1,9);
      book(_e["K3_45"],12,1,10);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //apply the final state of all particles and the one for only charged particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event,"CFS");
      const Particles& particles = cfs.particles();

      //because it does not make sense to have a collision with the numparticles is less than two,we use the vetoEvent so that if there is an event like this it does not run the analysis and goes to the next one

      if (numParticles<2) {
        MSG_DEBUG("Failed leptonic cut");
        vetoEvent;
      }

      //apply DIS kinematics in the event
      const DISKinematics& dk = apply<DISKinematics>(event,"Kinematics");
      //const DISLepton& dl = apply<DISLepton>(event,"Lepton");

      //get the DIS Kinematics but not in a loop because they are variables that descrube the type of event not the particles
      double Q2 = dk.Q2();
      double W2 = dk.W2();
      double W = std::sqrt(W2);
      bool cut1 = W<220. && W>80.;
      if (!cut1) vetoEvent;
      bool cut = Q2<1000. && Q2>10. && W>80. && W<220.;
      if (!cut) vetoEvent;

      double Efwd = 0. ;
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];
        double theta = p.theta()/degree ;
        if ( inRange(theta,4.4,15.) ) {
          Efwd = Efwd + p.E() ;
        }
      }

      bool cut_fwd = Efwd > 0.5 && dk.beamLepton().E() > 12. ;
      if (!cut_fwd) vetoEvent ;

      const size_t idx = _g["mult1"]->binAt(W).index();
      if (0 < idx && idx < 5) {
        _Nevt_after_cuts[idx-1]->fill();
      }


      //boost to the hadronic centre of mass frame
      const LorentzTransform hcmboost = dk.boostHCM();
      int kall = 0.;
      int k1 = 0.;
      int k2 = 0.;
      int k3 = 0.;
      int k4 = 0.;
      int k10 = 0.;
      int k11 = 0.;
      int k12 = 0.;
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        const double etahcm_char = hcmMom.eta();
        if (etahcm_char>0.)  ++kall;
        if (etahcm_char>1. && etahcm_char<2.) {
          ++k1;
        }
        if (etahcm_char>1. && etahcm_char<3.) {
          ++k2;
        }
        if (etahcm_char>1. && etahcm_char<4.) {
          ++k3;
        }
        if (etahcm_char>1. && etahcm_char<5.) {
          ++k4;
        }
        if (etahcm_char>2. && etahcm_char<3.) {
          ++k10;
        }
        if (etahcm_char>3. && etahcm_char<4.) {
          ++k11;
        }
        if (etahcm_char>4. && etahcm_char<5.) {
          ++k12;
        }

      }
      // cout<<k1<<endl;
      _g["mult_all"]->fill(W,kall);
      _g["mult10_all"]->fill(W,k10);
      _g["mult11_all"]->fill(W,k11);
      _g["mult12_all"]->fill(W,k12);
      _g["mult1"]->fill(W,k1);
      _g["mult2"]->fill(W,k2);
      _g["mult3"]->fill(W,k3);
      _g["mult4"]->fill(W,k4);

    }

    /// Normalise histograms etc., after the run
    void finalize() {
      int iW = 0 ;
      double iq  ;
      double mean, dispersion, cq, R2, R3, K3;
      for (auto& histo : _g["mult_all"]->bins()) {
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion, cq, R2, R3, K3);

        // just to have some values, needs to be corrected
        _e["mean0"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_0"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_0"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_0"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_0"]->bin(iW+1).set(R3, safediv(R3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_0"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_0"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_0"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_0"]->bin(iW+1).set(cq, safediv(cq, mean));

        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;
      for (auto& histo : _g["mult1"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);

        _e["mean12"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_12"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_12"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_12"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_12"]->bin(iW+1).set(R3, safediv(R3, mean));
        _e["K3_12"]->bin(iW+1).set(K3, safediv(K3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_12"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_12"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_12"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_12"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (auto& histo : _g["mult2"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq, mean, dispersion,cq, R2, R3, K3);

        _e["mean13"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_13"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_13"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_13"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_13"]->bin(iW+1).set(R3, safediv(R3, mean));
        _e["K3_13"]->bin(iW+1).set(K3, safediv(K3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_13"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_13"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_13"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_13"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (auto& histo : _g["mult3"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["mean14"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_14"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_14"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_14"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_14"]->bin(iW+1).set(R3, safediv(R3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_14"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_14"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3) ;
        _e["D4_14"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_14"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (auto& histo : _g["mult4"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["mean15"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_15"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_15"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_15"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_15"]->bin(iW+1).set(R3, safediv(R3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_15"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_15"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_15"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_15"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (auto& histo : _g["mult10_all"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["mean23"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_23"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_23"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_23"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_23"]->bin(iW+1).set(R3, safediv(R3, mean));
        _e["K3_23"]->bin(iW+1).set(K3, safediv(K3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_23"]->bin(iW+1).set(dispersion, dispersion/mean);
        _e["C3_23"]->bin(iW+1).set(cq, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_23"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_23"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }

      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (auto& histo : _g["mult11_all"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion, cq, R2, R3, K3);
        _e["mean34"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_34"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_34"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_34"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_34"]->bin(iW+1).set(R3, safediv(R3, mean));
        _e["K3_34"]->bin(iW+1).set(K3, safediv(K3, mean));
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_34"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_34"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_34"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_34"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (auto& histo : _g["mult12_all"]->bins()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        iq = 2 ;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["mean45"]->bin(iW+1).set(mean, 0.5*dispersion);
        _e["D2_45"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C2_45"]->bin(iW+1).set(cq, safediv(cq, mean));
        _e["R2_45"]->bin(iW+1).set(R2, safediv(R2, mean));
        _e["R3_45"]->bin(iW+1).set(R3, safediv(R3, mean));
        _e["K3_45"]->bin(iW+1).set(K3,K3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D3_45"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C3_45"]->bin(iW+1).set(cq, safediv(cq, mean));
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq, mean, dispersion, cq, R2, R3, K3);
        _e["D4_45"]->bin(iW+1).set(dispersion, safediv(dispersion, mean));
        _e["C4_45"]->bin(iW+1).set(cq, safediv(cq, mean));
        ++iW;
      }

    }


    inline void _histo_to_moments(const BinnedHistoPtr<int>& histo_input,
                                  double iq, double& mean, double& dispersion, double& cq,
                                  double& R2, double& R3, double& K3) {

      double mysumWX = 0.;
      double mysumW = 0.;

      if (histo_input->effNumEntries() == 0 || histo_input->sumW() == 0) {
        MSG_WARNING("Requested mean of a distribution with no net fill weights");
      }
      else {
        // loop to calcualte mean
        for (auto& b : histo_input->bins()) { // loop over points
          mysumWX  += b.sumW()*b.xEdge();
          mysumW   += b.sumW();
        }
        mean = safediv(mysumWX, mysumW);

        // loop to calculate dispersion (variance)
        double var = 0., c = 0., r2 = 0., r3 = 0.;
        for (auto& b : histo_input->bins()) { // loop over points
          const double xval = b.xEdge();
          const double weight = b.sumW();
          var = var + weight * pow((xval - mean),iq) ;
          c = c+pow(xval,iq)*weight;
          r2 = r2+xval*(xval-1)*weight;
          r3 = r3+xval*(xval-1)*(xval-2)*weight;
        }
        var = safediv(var, mysumW);
        dispersion = pow(var,1./iq);
        cq = safediv(c, mysumW*pow(mean,iq));
        R2 = safediv(r2, mysumW*pow(mean,2));
        R3 = safediv(r3, mysumW*pow(mean,3));
        K3 = R3 - 3*R2 + 2;
      }
    }

    ///@}


  private:

    /// @name Histograms
    ///@{

    map<string,HistoGroupPtr<double,int>> _g;
    map<string,Estimate1DPtr> _e;
    CounterPtr _Nevt_after_cuts[4];

    ///@}

  };

  RIVET_DECLARE_PLUGIN(H1_1996_I422230);

}
