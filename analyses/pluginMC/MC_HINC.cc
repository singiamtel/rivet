// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief MC validation analysis for higgs [-> tau tau] events
  class MC_HINC : public Analysis {
  public:

    /// Default constructor
    MC_HINC()
      : Analysis("MC_HINC")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETATAUMAX", 3.5);
      const double ptcut = getOption<double>("PTTAUMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;
      /// @todo Er, FS taus?
      DileptonFinder hfinder(125*GeV, 0.0, cut && Cuts::abspid == PID::TAU, Cuts::massIn(115*GeV, 135*GeV));

      declare(hfinder, "Hfinder");
      book(_h_H_mass ,"H_mass", 50, 119.7, 120.3);
      book(_h_H_pT ,"H_pT", logspace(100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));
      book(_h_H_pT_peak ,"H_pT_peak", 25, 0.0, 25.0);
      book(_h_H_y ,"H_y", 40, -4, 4);
      book(_h_H_phi ,"H_phi", 25, 0.0, TWOPI);
      book(_h_lepton_pT ,"lepton_pT", logspace(100, 10.0, 0.25*(sqrtS()>0.?sqrtS():14000.)/GeV));
      book(_h_lepton_eta ,"lepton_eta", 40, -4, 4);
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const DileptonFinder& hfinder = apply<DileptonFinder>(e, "Hfinder");
      if (hfinder.bosons().size() != 1) vetoEvent;

      FourMomentum hmom(hfinder.bosons()[0].momentum());
      _h_H_mass->fill(hmom.mass()/GeV);
      _h_H_pT->fill(hmom.pT()/GeV);
      _h_H_pT_peak->fill(hmom.pT()/GeV);
      _h_H_y->fill(hmom.rapidity());
      _h_H_phi->fill(hmom.phi());
      for (const Particle& l : hfinder.constituents()) {
        _h_lepton_pT->fill(l.pT()/GeV);
        _h_lepton_eta->fill(l.eta());
      }
    }


    /// Finalize
    void finalize() {
      const double xsec = crossSection()/picobarn;
      normalize(_h_H_mass, xsec);
      normalize(_h_H_pT, xsec);
      normalize(_h_H_pT_peak, xsec);
      normalize(_h_H_y, xsec);
      normalize(_h_H_phi, xsec);
      normalize(_h_lepton_pT, xsec);
      normalize(_h_lepton_eta, xsec);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_H_mass;
    Histo1DPtr _h_H_pT;
    Histo1DPtr _h_H_pT_peak;
    Histo1DPtr _h_H_y;
    Histo1DPtr _h_H_phi;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;
    /// @}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_HINC);

}
