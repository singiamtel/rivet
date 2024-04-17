// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Event-shape variables in deep-inelastic scattering at HERA
  class H1_2006_I699835 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2006_I699835);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs, "FS");

      const DISFinalState DISfs(DISFrame::BREIT);

      const FinalState DISfsCut(DISfs, Cuts::eta < 0);

      declare(Thrust(DISfsCut), "ThrustCut");

      declare(DISKinematics(), "Kinematics");

      //Book histograms for different Q ranges:
      book(_Nevt_after_cuts, "TMP/Nevt_after_cuts");

      const vector<double> QEdges{14., 16., 20., 30., 50., 70., 100., 200.};
      book(_h_tauc, QEdges);
      book(_h_tau, QEdges);
      book(_h_B, QEdges);
      book(_h_rho, QEdges);
      for (size_t iQ=0; iQ < QEdges.size()-1; ++iQ) {
        book(_h_tauc->bin(iQ+1), 1+iQ,1,1);
        book(_h_tau->bin(iQ+1),  8+iQ,1,1);
        book(_h_B->bin(iQ+1),   15+iQ,1,1);
        book(_h_rho->bin(iQ+1), 22+iQ,1,1);
        book(_Nevt_after_cuts_Q[iQ], "TMP/Nevt_after_cuts_Q"+ to_string(iQ));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");

      // The kinematic region covered by the analysis is defined by ranges of Q² and y
      if (dk.Q2() < 196 or dk.Q2() > 40000 or dk.y() < 0.1 or dk.y() > 0.7) {
      	vetoEvent;
      }

      const double Q = sqrt(dk.Q2());
      _Nevt_after_cuts->fill();
      const size_t iQ = _h_tauc->binAt(Q).index();
      if (0 < iQ && iQ < 8)  _Nevt_after_cuts_Q[iQ-1]->fill();


      // Boost to hadronic Breit frame
      const LorentzTransform breitboost = dk.boostBreit();

      const FinalState& fs = apply<FinalState>(event, "FS");

      /*Calculate event shape variables:
      thrust_num is \sum |pz_h|
      thrust_den is \sum |p_h|
      b_num is \sum |pt_h|
      sumE is the sum of energies
      (All sums run through the particles in the current hemisphere)
      */
      double thrust_num, thrust_den, b_num, sumE;
      thrust_num = thrust_den = b_num = sumE = 0;
      Vector3 sumMom;

      for (const Particle& p : fs.particles()) {
      	// Boost to Breit frame
        const FourMomentum breitMom = breitboost.transform(p.momentum());
      	if (breitMom.eta() < 0) {
          thrust_num += abs(breitMom.pz());
          thrust_den += breitMom.p();
          b_num += abs(breitMom.pt());
          sumMom.operator+=(breitMom.p3());
          sumE += breitMom.E();
      	}
      }

      //The energy in the current hemisphere must exceed a certain value.
      if (sumE <= Q/10.0)  vetoEvent;

      /* Comment from A. Galvan:
      Thrust here is with respect to the z axis (tau in the paper), while the
      one from Rivet projection is with respect to the maximum thrust
      axis (1 - tau_c in the paper)
      */

      double thrust_mine = 1.0 - ((double)thrust_num)/((double)thrust_den);
      double b_mine = ((double)b_num)/(2.0*(double)thrust_den);
      double rho_num = thrust_den*thrust_den - sumMom.dot(sumMom);
      double rho_mine = ((double)rho_num)/(4.0*(double)thrust_den*(double)thrust_den);

      const Thrust& thrCut = apply<Thrust>(event, "ThrustCut");

      //Fill histograms:
      _h_tauc->fill(Q, 1.0 - thrCut.thrust());
      _h_tau->fill(Q, thrust_mine);
      _h_B->fill(Q, b_mine);
      _h_rho->fill(Q, rho_mine);

      /* Comment from A. Galvan:
       As for the C-parameter, my results did not fit the reference data at
       all. The formula given in the paper is a bit ambiguous, because there is
       a sum that runs through all pairs of particles h,h' and I was not sure
       whether each pair should be counted twice (h,h' and h',h) or not, or if
       the pair of a particle with itself (hh) should be considered.
       That is why I tried all of the possibilities, but none of them worked.
      */
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // need to multiply by bin widths
      for (size_t iQ=0; iQ < _h_tauc->numBins(); ++iQ) {
        const double Nev = dbl(*_Nevt_after_cuts_Q[iQ]);
        if (Nev != 0) {
          scale(_h_tauc->bin(iQ+1), 1./Nev);
          for(auto & bin : _h_tauc->bin(iQ+1)->bins())
            bin.scaleW(bin.xWidth());
          scale(_h_tau->bin(iQ+1), 1./Nev);
          for(auto & bin : _h_tau->bin(iQ+1)->bins())
            bin.scaleW(bin.xWidth());
          scale(_h_B->bin(iQ+1), 1./Nev);
          for(auto & bin : _h_B->bin(iQ+1)->bins())
            bin.scaleW(bin.xWidth());
          scale(_h_rho->bin(iQ+1), 1./Nev);
          for(auto & bin : _h_rho->bin(iQ+1)->bins())
            bin.scaleW(bin.xWidth());
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DGroupPtr _h_tauc, _h_tau, _h_B, _h_rho;
    CounterPtr _Nevt_after_cuts;
    CounterPtr _Nevt_after_cuts_Q[7];
    ///@}

  };


  RIVET_DECLARE_PLUGIN(H1_2006_I699835);

}
