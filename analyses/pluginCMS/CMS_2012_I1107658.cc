// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Underlying event activity in the Drell-Yan process at 7 TeV
  class CMS_2012_I1107658 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I1107658);


    /// Initialization
    void init() {

      /// @note Using a bare muon Z (but with a clustering radius!?)
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      DileptonFinder zfinder(91.2*GeV, 0.2, cut && Cuts::abspid == PID::MUON, Cuts::massIn(4*GeV, 140*GeV));
      declare(zfinder, "DileptonFinder");

      ChargedFinalState nonmuons(Cuts::abseta < 2 && Cuts::pT > 500*MeV && Cuts::abspid != PID::MUON);
      declare(nonmuons, "nonmuons");

      book(_h_Nchg_towards_pTmumu                 ,1, 1, 1);
      book(_h_Nchg_transverse_pTmumu              ,2, 1, 1);
      book(_h_Nchg_away_pTmumu                    ,3, 1, 1);
      book(_h_pTsum_towards_pTmumu                ,4, 1, 1);
      book(_h_pTsum_transverse_pTmumu             ,5, 1, 1);
      book(_h_pTsum_away_pTmumu                   ,6, 1, 1);
      book(_h_avgpT_towards_pTmumu                ,7, 1, 1);
      book(_h_avgpT_transverse_pTmumu             ,8, 1, 1);
      book(_h_avgpT_away_pTmumu                   ,9, 1, 1);
      book(_h_Nchg_towards_plus_transverse_Mmumu  ,10, 1, 1);
      book(_h_pTsum_towards_plus_transverse_Mmumu ,11, 1, 1);
      book(_h_avgpT_towards_plus_transverse_Mmumu ,12, 1, 1);
      book(_h_Nchg_towards_zmass_81_101           ,13, 1, 1);
      book(_h_Nchg_transverse_zmass_81_101        ,14, 1, 1);
      book(_h_Nchg_away_zmass_81_101              ,15, 1, 1);
      book(_h_pT_towards_zmass_81_101             ,16, 1, 1);
      book(_h_pT_transverse_zmass_81_101          ,17, 1, 1);
      book(_h_pT_away_zmass_81_101                ,18, 1, 1);
      book(_h_Nchg_transverse_zpt_5               ,19, 1, 1);
      book(_h_pT_transverse_zpt_5                 ,20, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const DileptonFinder& zfinder = apply<DileptonFinder>(event, "DileptonFinder");

      if (zfinder.bosons().size() != 1) vetoEvent;

      double Zpt = zfinder.bosons()[0].pT()/GeV;
      double Zphi = zfinder.bosons()[0].phi();
      double Zmass = zfinder.bosons()[0].mass()/GeV;

      Particles particles = apply<VetoedFinalState>(event, "nonmuons").particles();

      int nTowards = 0;
      int nTransverse = 0;
      int nAway = 0;
      double ptSumTowards = 0;
      double ptSumTransverse = 0;
      double ptSumAway = 0;

      for (const Particle& p : particles) {
        double dphi = fabs(deltaPhi(Zphi, p.phi()));
        double pT = p.pT();

        if ( dphi < M_PI/3 ) {
          nTowards++;
          ptSumTowards += pT;
          if (Zmass > 81. && Zmass < 101.) _h_pT_towards_zmass_81_101->fill(pT);
        } else if ( dphi < 2.*M_PI/3 ) {
          nTransverse++;
          ptSumTransverse += pT;
          if (Zmass > 81. && Zmass < 101.) _h_pT_transverse_zmass_81_101->fill(pT);
          if (Zpt < 5.) _h_pT_transverse_zpt_5->fill(pT);
        } else {
          nAway++;
          ptSumAway += pT;
          if (Zmass > 81. && Zmass < 101.) _h_pT_away_zmass_81_101->fill(pT);
        }

      } // Loop over particles


      const double area = 8./3.*M_PI;
      if (Zmass > 81. && Zmass < 101.) {
        _h_Nchg_towards_pTmumu->         fill(Zpt, 1./area * nTowards);
        _h_Nchg_transverse_pTmumu->      fill(Zpt, 1./area * nTransverse);
        _h_Nchg_away_pTmumu->            fill(Zpt, 1./area * nAway);
        _h_pTsum_towards_pTmumu->        fill(Zpt, 1./area * ptSumTowards);
        _h_pTsum_transverse_pTmumu->     fill(Zpt, 1./area * ptSumTransverse);
        _h_pTsum_away_pTmumu->           fill(Zpt, 1./area * ptSumAway);
        if (nTowards > 0)    _h_avgpT_towards_pTmumu->    fill(Zpt, ptSumTowards/nTowards);
        if (nTransverse > 0) _h_avgpT_transverse_pTmumu-> fill(Zpt, ptSumTransverse/nTransverse);
        if (nAway > 0)       _h_avgpT_away_pTmumu->       fill(Zpt, ptSumAway/nAway);
        _h_Nchg_towards_zmass_81_101->   fill(nTowards);
        _h_Nchg_transverse_zmass_81_101->fill(nTransverse);
        _h_Nchg_away_zmass_81_101->      fill(nAway);
      }

      if (Zpt < 5.) {
        _h_Nchg_towards_plus_transverse_Mmumu->fill(Zmass, (nTowards + nTransverse)/(2.*area));
        _h_pTsum_towards_plus_transverse_Mmumu->fill(Zmass, (ptSumTowards + ptSumTransverse)/(2.*area));
        if ((nTowards + nTransverse) > 0) _h_avgpT_towards_plus_transverse_Mmumu->fill(Zmass, (ptSumTowards + ptSumTransverse)/(nTowards + nTransverse));
        _h_Nchg_transverse_zpt_5->fill(nTransverse);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pT_towards_zmass_81_101,    safediv(1, _h_Nchg_towards_zmass_81_101->integral(), 0));
      scale(_h_pT_transverse_zmass_81_101, safediv(1, _h_Nchg_transverse_zmass_81_101->integral(), 0));
      scale(_h_pT_away_zmass_81_101,       safediv(1, _h_Nchg_away_zmass_81_101->integral(), 0));
      scale(_h_pT_transverse_zpt_5,        safediv(1, _h_Nchg_transverse_zpt_5->integral(), 0));
      normalize(_h_Nchg_towards_zmass_81_101);
      normalize(_h_Nchg_transverse_zmass_81_101);
      normalize(_h_Nchg_away_zmass_81_101);
      normalize(_h_Nchg_transverse_zpt_5);
    }


  private:


    /// @name Histogram objects
    /// @{
    Profile1DPtr _h_Nchg_towards_pTmumu;
    Profile1DPtr _h_Nchg_transverse_pTmumu;
    Profile1DPtr _h_Nchg_away_pTmumu;
    Profile1DPtr _h_pTsum_towards_pTmumu;
    Profile1DPtr _h_pTsum_transverse_pTmumu;
    Profile1DPtr _h_pTsum_away_pTmumu;
    Profile1DPtr _h_avgpT_towards_pTmumu;
    Profile1DPtr _h_avgpT_transverse_pTmumu;
    Profile1DPtr _h_avgpT_away_pTmumu;
    Profile1DPtr _h_Nchg_towards_plus_transverse_Mmumu;
    Profile1DPtr _h_pTsum_towards_plus_transverse_Mmumu;
    Profile1DPtr _h_avgpT_towards_plus_transverse_Mmumu;
    Histo1DPtr _h_Nchg_towards_zmass_81_101;
    Histo1DPtr _h_Nchg_transverse_zmass_81_101;
    Histo1DPtr _h_Nchg_away_zmass_81_101;
    Histo1DPtr _h_pT_towards_zmass_81_101;
    Histo1DPtr _h_pT_transverse_zmass_81_101;
    Histo1DPtr _h_pT_away_zmass_81_101;
    Histo1DPtr _h_Nchg_transverse_zpt_5;
    Histo1DPtr _h_pT_transverse_zpt_5;
    /// @}

  };


  // Hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2012_I1107658);

}
