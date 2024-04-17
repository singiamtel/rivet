// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Transverse momentum spectra of charged particles in DIS
  class H1_1996_I424463 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1996_I424463);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Book projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(ChargedFinalState(), "CFS");
      declare(FinalState(), "FS");

      // Book histograms
      book(_NevAll, "TMP/Nev_all");
      int ixx = 0 ;
      for (size_t ix = 0; ix < 10; ++ix) {
        book(_Nevt_after_cuts[ix], "TMP/Nevt_after_cuts" + to_string(ix));
        for(unsigned int ih=0;ih<2;++ih) {
          book(_h_dndpt_eta_bin[ih][ix], ih*10 + ix+1 , 1, 1);
          book(_hdndeta_bin    [ih][ix], ih*10 + ix+29, 1, 1);
        }
        if (ix != 6 && ix != 9) {
          book(_hdndptmax_low_eta_bin[ixx], ixx+21, 1, 1);
          ixx=ixx+1;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edgespT.empty()) {
        _edgespT    = _h_dndpt_eta_bin   [0][0]->xEdges();
        _edgesEta   = _hdndeta_bin       [0][0]->xEdges();
        _edgespTMax = _hdndptmax_low_eta_bin[0]->xEdges();
      }
      const FinalState& fs = apply<FinalState>(event, "FS");
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;
      double W2 = dk.W2()/GeV;

      // Momentum of the scattered lepton
      FourMomentum leptonMom = dl.out().momentum();
      double enel = leptonMom.E()/GeV;
      double thel = 180.-leptonMom.angle(dl.in().momentum())/degree;

      _NevAll -> fill() ;

      bool cut = y > 0.05 && Q2 > 5. && Q2 < 100.&& enel > 12. && W2 > 4400. && thel > 157. && thel < 173.;
      if (!cut) vetoEvent;

      int ibin[10] ;
      for (int i = 0; i < 10; i++) ibin[i] = 0;
      if ( 5. < Q2 && Q2 < 50. && x > 0.0001 && x < 0.0010) ibin[0] = 1;
      if ( 5. < Q2 && Q2 < 10. && x > 0.0001 && x < 0.0002) ibin[1] = 1;
      if ( 6. < Q2 && Q2 < 10. && x > 0.0002 && x < 0.0005) ibin[2] = 1;
      if (10. < Q2 && Q2 < 20. && x > 0.0002 && x < 0.0005) ibin[3] = 1;
      if (10. < Q2 && Q2 < 20. && x > 0.0005 && x < 0.0008) ibin[4] = 1;
      if (10. < Q2 && Q2 < 20. && x > 0.0008 && x < 0.0015) ibin[5] = 1;
      if (10. < Q2 && Q2 < 20. && x > 0.0015 && x < 0.0040) ibin[6] = 1;
      if (20. < Q2 && Q2 < 50. && x > 0.0005 && x < 0.0014) ibin[7] = 1;
      if (20. < Q2 && Q2 < 50. && x > 0.0014 && x < 0.0030) ibin[8] = 1;
      if (20. < Q2 && Q2 < 50. && x > 0.0030 && x < 0.0100) ibin[9] = 1;
      for (int i = 0; i < 10; i++) {
        if (ibin[i]==1) _Nevt_after_cuts[i]->fill();
      }

      // Extract the particles other than the lepton
      /// @todo Improve to avoid HepMC digging
      Particles particles;
      particles.reserve(fs.particles().size());
      ConstGenParticlePtr dislepGP = dl.out().genParticle();
      for (const Particle& p : fs.particles()) {
        ConstGenParticlePtr loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }

      // Boost to hadronic CM
      const LorentzTransform hcmboost = dk.boostHCM();

      // Loop over the particles
      int mult = 0 ;
      double ptmax_high[10], ptmax_low[10];
      for (int i = 0; i < 10; i++) {
        ptmax_high[i] = 0.; ptmax_low[i] = 0.;
      }
      double EtSum = 0;
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];

        double eta = p.momentum().pseudorapidity();
        // Boost to hcm
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());

        // Apply safety cuts
        if (eta > -5  && eta < 10.) {
          mult = mult + 1;
          double pThcm =hcmMom.pT() ;
          double etahcm = hcmMom.pseudorapidity();
          if (etahcm > 0. && etahcm < 2.0) {
            EtSum = EtSum + hcmMom.Et();
          }
          if (PID::charge(p.pid()) != 0) {
            if (etahcm > 0.5 && etahcm < 1.5) {
              for ( int i=0; i< 10; i++) {
                if (ibin[i]==1) {
                  fillpT(1,i,pThcm);
                  if (pThcm > ptmax_low[i] ) ptmax_low[i] = pThcm;
                }
              }

            }
            if (etahcm > 1.5 && etahcm < 2.5){
              for ( int i=0; i< 10; i++) {
                if (ibin[i]==1) {
                  fillpT(0,i,pThcm);
                  if (pThcm > ptmax_high[i] ) ptmax_high[i] = pThcm;
                }
              }
            }
            for ( int i=0; i< 10; i++) {
              if (ibin[i]==1)               fillEta(1,i,etahcm);
              if (ibin[i]==1 && pThcm > 1.) fillEta(0,i,etahcm);
            }
          }
        }  // end of loop over the particles
      }
      int ii=0;
      for ( int i=0; i< 10; i++) {
        if (i != 6 && i != 9 ) {
          if ( ibin[i]==1 && EtSum > 6. ) {
            fillpTMax(ii,ptmax_low[i]);
          }
          ii=ii+1;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_DEBUG("All events: " << _NevAll->val() << " after cuts: "<<  _Nevt_after_cuts[0]->val());
      MSG_DEBUG("Cut1 events: " << _NevAll->val() << " after cuts: "<<  _Nevt_after_cuts[1]->val());
      int ii = 0;
      for ( int i=0; i< 10; i++) {
        if (_Nevt_after_cuts[i]->val()  != 0) {
          for(unsigned int ih=0;ih<2;++ih) {
            scale(_h_dndpt_eta_bin[ih][i], 1./ *_Nevt_after_cuts[i]);
            size_t ioff=0;
            if( (ih==0 && (i==6 || i==9)) ||
                (ih==1 && (i<=4 || i==7)))    ioff=1;
            else if (ih==1 && (i==5 || i==8)) ioff=2;
            for(auto & b : _h_dndpt_eta_bin[ih][i]->bins()) {
              const size_t idx = b.index()+ioff;
              b.scaleW(1./_axispT.width(idx));
            }
            scale(_hdndeta_bin    [ih][i], 1./ *_Nevt_after_cuts[i]);
            ioff=0;
            if( (ih==0 && (i==5 || i==6 || i==8 || i==9)) ||
                (ih==1 && (i==0 || i==4 || i==5 || i==8))) ioff=1;
            else if (ih==1 && (i==6 || i==9))       ioff=2;
            for(auto & b : _hdndeta_bin[ih][i]->bins()) {
              const size_t idx = b.index()+ioff;
              b.scaleW(1./_axisEta.width(idx));
            }
          }
        }
        if ( i !=6 && i!=9) {
          if (_Nevt_after_cuts[i]->val()  != 0) {
            normalize(_hdndptmax_low_eta_bin[ii]);
            for(auto & b : _hdndptmax_low_eta_bin[ii]->bins()) {
              const size_t idx = b.index();
              b.scaleW(1./_axispTMax.width(idx));
            }
          }
          ii=ii+1;
        }
      }
      
    }
    /// @}

    void fillEta(const unsigned int ix, const unsigned int iy, const double value) {
      string edge = "OTHER";
      const size_t idx = _axisEta.index(value);
      if (idx && idx <= _edgesEta.size()) {
        if ( ( (ix==0 && iy>=1&& iy<=3) || (ix==0 && (iy==1 || iy==3)) ) && idx == _edgesEta.size())
          ;
        else if (idx==1 && ( (ix==0 && (iy==5 || iy==6 || iy==8 || iy==9)) ||
                             (ix==1 && (iy==0 || iy==4 || iy==5 || iy==8))))
          ;
        else if (idx<=2 && (ix==1 && (iy==6 || iy==9)) )
          ;
        else {
          edge = _edgesEta[idx-1];
        }
      }
      _hdndeta_bin[ix][iy]->fill(edge);
    }
    
    void fillpT(const unsigned int ix, const unsigned int iy, const double value) {
      string edge = "OTHER";
      const size_t idx = _axispT.index(value);
      if (idx && idx <= _edgespT.size()) {
        if(((ix==0 && iy==8) || (ix==1 && (iy==2 || iy==4 || iy==6))) && idx == _edgespT.size())
          ;
        else if (idx==1 && ( (ix==0 && (iy==6 || iy==9)) ||
                             (ix==1 && (iy<=4 || iy==7))))
          ;
        else if (idx<=2 && (ix==1 && (iy==5 || iy==8)))
          ;
        else
          edge = _edgespT[idx-1];
      }
      _h_dndpt_eta_bin[ix][iy]->fill(edge);
    }

    void fillpTMax(const unsigned int ix, const double value) {
      string edge = "OTHER";
      const size_t idx = _axispTMax.index(value);
      if (idx && idx <= _edgespTMax.size()) {
        if(idx==_edgespTMax.size() && (ix==1 || ix==2 || ix==4))
          ;
        else if(idx>=_edgespTMax.size()-1 && ix==5)
          ;
        else
          edge = _edgespTMax[idx-1];
      }
      _hdndptmax_low_eta_bin[ix]->fill(edge);
    }

    /// @name Histograms
    /// @{
    array<CounterPtr,10> _Nevt_after_cuts;
    BinnedHistoPtr<string> _h_dndpt_eta_bin[2][10], _hdndeta_bin[2][10], _hdndptmax_low_eta_bin[8];
    CounterPtr _NevAll ;

    vector<string> _edgespT,_edgesEta,_edgespTMax;
    YODA::Axis<double> _axispT    = YODA::Axis<double>{ 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.825, 2.125, 2.525, 3.125, 4.0, 5.0};
    YODA::Axis<double> _axisEta   = YODA::Axis<double>{ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.00};
    YODA::Axis<double> _axispTMax = YODA::Axis<double>{0.005, 0.255, 0.505, 0.755, 1.005, 1.255, 1.505, 1.755, 2.065, 2.5, 3.125, 4.0, 5.0};
    /// @}

  };


  RIVET_DECLARE_ALIASED_PLUGIN(H1_1996_I424463, H1_1997_I424463);

}
