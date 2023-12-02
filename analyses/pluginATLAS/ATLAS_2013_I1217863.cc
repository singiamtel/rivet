// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Electroweak Wjj production at 8 TeV
  class ATLAS_2013_I1217863 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1217863);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get modes from the option system
      _mode = 2;
      _doZ  = true;
      _doW  = true;
      if ( getOption("LMODE") == "EL" ) { _mode = 2;}
      if ( getOption("LMODE") == "MU" ) _mode = 3;
      if ( getOption("LMODE") == "ZEL" ) {
        _mode = 2;
        _doW  = false;
      }
      if ( getOption("LMODE") == "ZMU" ) {
        _mode = 3;
        _doW  = false;
      }
      if ( getOption("LMODE") == "WEL" ) {
        _mode = 2;
        _doZ  = false;
      }
      if ( getOption("LMODE") == "WMU" ) {
        _mode = 3;
        _doZ  = false;
      }

      Cut cuts = Cuts::abseta < 2.47 && Cuts::pT > 25*GeV;
      VetoedFinalState jet_fs;

      // Z finder
      if (_doZ) {
        DileptonFinder zf(91.2*GeV, 0.1, cuts && Cuts::abspid == (_mode == 3 ? PID::MUON : PID::ELECTRON), Cuts::massIn(40.0*GeV, 1000.0*GeV));
        declare(zf, "ZF");
        jet_fs.addVetoOnThisFinalState(zf);
      }

      if (_doW) {
        // W finder
        declare("MET", MissingMomentum());
        LeptonFinder lf(0.1, cuts && Cuts::abspid == (_mode==3? PID::MUON : PID::ELECTRON));
        declare(lf, "Leptons");
        jet_fs.addVetoOnThisFinalState(lf);
      }

      // Leading photon
      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 2.37 && Cuts::pT > 15*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");
      jet_fs.addVetoOnThisFinalState(photonfs);

      // Jets
      FastJets jets(jet_fs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::NONE);
      declare(jets, "Jets");

      // FS excluding only the leading photon
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "isolatedFS");


      // Book histograms
      if (_doZ) {
        book(_hist_EgammaT_inclZ   ,11, 1, _mode); // dSigma / dE^gamma_T for Njet >= 0
        book(_hist_EgammaT_exclZ   ,12, 1, _mode); // dSigma / dE^gamma_T for Njet = 0
        book(_hist_Njet_EgammaT15Z ,17, 1, _mode); // dSigma / dNjet for E^gamma_T >= 15
        book(_hist_Njet_EgammaT60Z ,18, 1, _mode); // dSigma / dNjet for E^gamma_T >= 60
        book(_hist_mZgamma         ,20, 1, _mode); // dSigma / dm^{Zgamma}
      }
      if (_doW){
        book(_hist_EgammaT_inclW   , 7, 1, _mode); // dSigma / dE^gamma_T for Njet >= 0
        book(_hist_EgammaT_exclW   , 8, 1, _mode); // dSigma / dE^gamma_T for Njet = 0
        book(_hist_Njet_EgammaT15W ,15, 1, _mode); // dSigma / dNjet for E^gamma_T >= 15
        book(_hist_Njet_EgammaT60W ,16, 1, _mode); // dSigma / dNjet for E^gamma_T >= 60
        book(_hist_mWgammaT        ,19, 1, _mode); // dSigma / dm^{Zgamma}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve leading photon
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1)  vetoEvent;
      const Particle& leadingPhoton = photons[0];
      if (leadingPhoton.Et() < 15.0*GeV) vetoEvent;
      if (leadingPhoton.abseta() > 2.37) vetoEvent;

      // Check photon isolation
      double coneEnergy(0.0);
      Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
      for (const Particle& p : fs) {
        if ( deltaR(leadingPhoton, p) < 0.4 )  coneEnergy += p.E();
      }
      if (coneEnergy / leadingPhoton.E() >= 0.5 )  vetoEvent;

      if (_doW) {
	// Retrieve W boson candidate
	const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
	if (pmiss.pT() > 35*GeV) {

	  const Particles& ls = apply<LeptonFinder>(event, "Leptons").particles();
	  const int ifound = closestMatchIndex(ls, pmiss, Kin::mass, 80.4*GeV);
	  if (ifound >= 0) {

	    // Retrieve constituent lepton
	    const Particle& lepton = ls[ifound];
	    if ( lepton.pT() > 25.0*GeV && lepton.abseta() < 2.47 ) { //< redundant cut

	      // Check photon-lepton overlap
	      if ( deltaR(leadingPhoton, lepton) > 0.7 ) {

		// Count jets
		const FastJets& jetfs = apply<FastJets>(event, "Jets");
		Jets jets = jetfs.jets(cmpMomByEt);
		int goodJets = 0;
		for (const Jet& j : jets) {
		  if ( !(j.Et() > 30.0*GeV) )  break;
		  if ( (j.abseta() < 4.4) &&				\
		       (deltaR(leadingPhoton, j) > 0.3) &&		\
		       (deltaR(lepton,        j) > 0.3) )  ++goodJets;
		}

		double Njets = double(goodJets) + 0.5;
		double photonEt = leadingPhoton.Et()*GeV;

		const FourMomentum& lep_gamma = lepton.momentum() + leadingPhoton.momentum();
		double term1 = sqrt(lep_gamma.mass2() + lep_gamma.pT2()) + pmiss.Et();
		double term2 = (lep_gamma + pmiss).pT2();
		double mWgammaT = sqrt(term1 * term1 - term2) * GeV;

		_hist_EgammaT_inclW->fill(photonEt);

		_hist_Njet_EgammaT15W->fill(Njets);

		if ( !goodJets )  _hist_EgammaT_exclW->fill(photonEt);

		if (photonEt > 40.0*GeV) {
		  _hist_mWgammaT->fill(mWgammaT);
		  if (photonEt > 60.0*GeV)  _hist_Njet_EgammaT60W->fill(Njets);
		}
	      }
	    }
	  }
	}
      }

      if (_doZ ){

	// Retrieve Z boson candidate
	const DileptonFinder& zf = apply<DileptonFinder>(event, "ZF");
	if ( zf.bosons().size() == 1 ) {
	  const Particle& Zboson  = zf.boson();
	  if ( (Zboson.mass() > 40.0*GeV) ) {

	    // Check charge of constituent leptons
	    const Particles& leptons = zf.constituents();
	    if (leptons.size() == 2 && leptons[0].charge() * leptons[1].charge() < 0.) {

	      bool lpass = true;
	      // Check photon-lepton overlap
	      for (const Particle& p : leptons) {
		if ( !(p.pT() > 25.0*GeV && p.abseta() < 2.47 && deltaR(leadingPhoton, p) > 0.7) )  lpass = false;
	      }
	      if ( lpass ) {

		// Count jets
		const FastJets& jetfs = apply<FastJets>(event, "Jets");
		Jets jets = jetfs.jets(cmpMomByEt);
		int goodJets = 0;
		for (const Jet& j : jets) {
		  if ( !(j.Et() > 30.0*GeV) )  break;
		  if ( (j.abseta() < 4.4) &&		    \
		       (deltaR(leadingPhoton, j) > 0.3) &&  \
		       (deltaR(leptons[0],    j) > 0.3) &&		\
		       (deltaR(leptons[1],    j) > 0.3) )  ++goodJets;
		}

		double Njets = double(goodJets) + 0.5;
		double photonEt = leadingPhoton.Et()*GeV;
		double mZgamma = (Zboson.momentum() + leadingPhoton.momentum()).mass() * GeV;

		_hist_EgammaT_inclZ->fill(photonEt);

		_hist_Njet_EgammaT15Z->fill(Njets);

		if ( !goodJets )   _hist_EgammaT_exclZ->fill(photonEt);

		if (photonEt >= 40.0*GeV) {
		  _hist_mZgamma->fill(mZgamma);
		  if (photonEt >= 60.0*GeV)  _hist_Njet_EgammaT60Z->fill(Njets);
		}
	      }
	    }
	  }
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs_fb = crossSection()/femtobarn;
      const double sumw = sumOfWeights();
      const double sf = xs_fb / sumw;

      if (_doZ) {
	scale(_hist_EgammaT_exclZ, sf);
	scale(_hist_EgammaT_inclZ, sf);
	normalize(_hist_Njet_EgammaT15Z);
	normalize(_hist_Njet_EgammaT60Z);
	normalize(_hist_mZgamma);
      }

      if (_doW) {
	scale(_hist_EgammaT_exclW, sf);
	scale(_hist_EgammaT_inclW, sf);
	normalize(_hist_Njet_EgammaT15W);
	normalize(_hist_Njet_EgammaT60W);
	normalize(_hist_mWgammaT);
      }

    }

    /// @}

  protected:

    size_t _mode;
    bool _doW;
    bool _doZ;

  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _hist_EgammaT_inclZ;
    Histo1DPtr _hist_EgammaT_exclZ;
    Histo1DPtr _hist_Njet_EgammaT15Z;
    Histo1DPtr _hist_Njet_EgammaT60Z;
    Histo1DPtr _hist_mZgamma;
    //
    Histo1DPtr _hist_EgammaT_inclW;
    Histo1DPtr _hist_EgammaT_exclW;
    Histo1DPtr _hist_Njet_EgammaT15W;
    Histo1DPtr _hist_Njet_EgammaT60W;
    Histo1DPtr _hist_mWgammaT;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1217863);

}
