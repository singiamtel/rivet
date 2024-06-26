// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// LHCb prompt charm hadron pT and rapidity spectra
  class LHCB_2016_I1490663 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1490663);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      declare(UnstableParticles(), "UFS");

      /// Book histograms
      book(_h_pdg411_Dplus_pT_y,     {2., 2.5, 3., 3.5, 4., 4.5});
      book(_h_pdg421_Dzero_pT_y,     {2., 2.5, 3., 3.5, 4., 4.5});
      book(_h_pdg431_Dsplus_pT_y,    {2., 2.5, 3., 3.5, 4., 4.5});
      book(_h_pdg413_Dstarplus_pT_y, {2., 2.5, 3., 3.5, 4., 4.5});
      for (size_t i = 1; i < _h_pdg411_Dplus_pT_y->numBins()+1; ++i) {
        size_t y = _h_pdg411_Dplus_pT_y->bin(i).index();
        book(_h_pdg411_Dplus_pT_y->bin(i),     1, 1, y);
        book(_h_pdg421_Dzero_pT_y->bin(i),     2, 1, y);
        book(_h_pdg431_Dsplus_pT_y->bin(i),    3, 1, y);
        book(_h_pdg413_Dstarplus_pT_y->bin(i), 4, 1, y);
      }

      book(_hbr_Dzero, {2., 2.5, 3., 3.5, 4., 4.5});
      book(_hbr_Dplus, {2., 2.5, 3., 3.5, 4., 4.5});
      book(_hbr_Ds,    {2., 2.5, 3., 3.5, 4., 4.5});
      book(_hbr_Dstar, {2., 2.5, 3., 3.5, 4., 4.5});
      for (size_t i = 1; i < _hbr_Dzero->numBins()+1; ++i) {
      	book(_hbr_Dzero->bin(i), "TMP/Dzero_b"+to_str(i), refData(9, 1, 2));
      	book(_hbr_Dplus->bin(i), "TMP/Dplus_b"+to_str(i), refData(9, 1, 2));
      	book(_hbr_Ds->bin(i),    "TMP/Ds_b"+to_str(i),    refData(9, 1, 2));
      	book(_hbr_Dstar->bin(i), "TMP/Dstar_b"+to_str(i), refData(9, 1, 2));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Use PrimaryHadrons to avoid double counting and automatically remove the contributions from unstable?
      const UnstableParticles &ufs = apply<UnstableParticles> (event, "UFS");
      for (const Particle& p : ufs.particles() ) {

        // We're only interested in charm hadrons
        //if (!p.isHadron() || !p.hasCharm()) continue;

        PdgId apid = p.abspid();

        // do not use Cuts::abspid to avoid supplemental iteration on particles?
        if ((apid != 411) && (apid != 421) && (apid != 431) && (apid != 413)) continue;

        // Experimental selection removes non-prompt charm hadrons: we ignore those from b decays
        if (p.fromBottom()) continue;

        // Kinematic acceptance
        const double y = p.absrap(); ///< Double analysis efficiency with a "two-sided LHCb"
        const double pT = p.pT()/GeV;

        // Fiducial acceptance of the measurements
        if ((pT > 10.0) || (y < 2.0) || (y > 4.5)) continue;

        Particles daus;

        switch (apid) {
        case 411:
          _h_pdg411_Dplus_pT_y->fill(y, pT);
          // veto on decay channel [D+ -> K- pi+ pi+]cc
          if (p.children().size() != 3) break;
          if ( ((p.children(Cuts::pid == -321).size() == 1) && (p.children(Cuts::pid == 211).size() == 2)) ||
          		 ((p.children(Cuts::pid == 321).size() == 1) && (p.children(Cuts::pid == -211).size() == 2)) )
          	_hbr_Dplus->fill(y, pT); // MSG_INFO("Found [ D+ -> K- pi+ pi+ ]cc..."); };
          break;
        case 421:
          _h_pdg421_Dzero_pT_y->fill(y, pT);
          // veto on decay channel [D0 -> K- pi+]cc
          if (p.children().size() != 2) break;
          if ( ((p.children(Cuts::pid == -321).size() == 1) && (p.children(Cuts::pid == 211).size() == 1)) ||
          		 ((p.children(Cuts::pid == 321).size() == 1) && (p.children(Cuts::pid == -211).size() == 1)) )
          	_hbr_Dzero->fill(y, pT); // MSG_INFO("Found [ D0 -> K- pi+ ]cc..."); };
          break;
        case 431:
          _h_pdg431_Dsplus_pT_y->fill(y, pT);
          //veto on decay channel [Ds+ -> [K+ K-]phi0 pi+]cc
          if (p.children().size() != 2) break;
          daus = p.children(Cuts::pid == 333);
          if ( (daus.size() == 1) && (p.children(Cuts::abspid == 211).size() == 1) &&
          		 (daus.front().children(Cuts::abspid ==321).size() == 2) )
          	_hbr_Ds->fill(y, pT); // MSG_INFO("Found [ Ds+ -> phi0(-> K+ K-) pi+ ]cc..."); };
          break;
        case 413:
          _h_pdg413_Dstarplus_pT_y->fill(y, pT);
          // veto on decay channel [D*+ -> [K- pi+]D0 pi+]cc
          if (p.children().size() != 2) break;
          daus = p.children(Cuts::pid == 421);
          if ( (daus.size() == 1) && (p.children(Cuts::abspid == 211).size() == 1) &&
          		( daus.front().children().size() == 2 ) &&
          		( ( (daus.front().children(Cuts::pid == -321).size() == 1 ) && (daus.front().children(Cuts::pid == 211).size() == 1 )	) ||
          		  ( (daus.front().children(Cuts::pid == 321).size() == 1 ) && (daus.front().children(Cuts::pid == -211).size() == 1 ) ) ) )
          	_hbr_Dstar->fill(y, pT); // MSG_INFO("Found [ D*+ -> D0 (-> K- pi+)cc pi+ ]cc..."); };
          break;
        default:
        	break;
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// Factor of 0.5 to correct for the abs(rapidity) used above
      const double scale_factor = 0.5 * crossSection()/microbarn / sumOfWeights();

      scale(_h_pdg411_Dplus_pT_y, scale_factor);
      scale(_h_pdg421_Dzero_pT_y, scale_factor);
      scale(_h_pdg431_Dsplus_pT_y, scale_factor);
      scale(_h_pdg413_Dstarplus_pT_y, scale_factor);

      // Do ratios
      for (int i = 0; i < 5; ++i) {
      	book(hr_DplusDzero[i], 9, 1, i+1);
      	book(hr_DsDzero[i], 10, 1, i+1);
      	book(hr_DstarDzero[i], 11, 1, i+1);
      	book(hr_DsDplus[i], 12, 1, i+1);
      	book(hr_DstarDplus[i], 13, 1, i+1);
      	book(hr_DsDstar[i], 14, 1, i+1);
      	divide(_hbr_Dplus->bin(i+1), _hbr_Dzero->bin(i+1), hr_DplusDzero[i]);
      	divide(_hbr_Ds->bin(i+1),    _hbr_Dzero->bin(i+1), hr_DsDzero[i]);
      	divide(_hbr_Dstar->bin(i+1), _hbr_Dzero->bin(i+1), hr_DstarDzero[i]);
      	divide(_hbr_Ds->bin(i+1),    _hbr_Dplus->bin(i+1), hr_DsDplus[i]);
      	divide(_hbr_Dstar->bin(i+1), _hbr_Dplus->bin(i+1), hr_DstarDplus[i]);
      	divide(_hbr_Ds->bin(i+1),    _hbr_Dstar->bin(i+1), hr_DsDstar[i]);
      	// scale 100x as measurement is in %
      	hr_DplusDzero[i]->scale(100.);
      	hr_DsDzero[i]->scale(100.);
      	hr_DstarDzero[i]->scale(100.);
      	hr_DsDplus[i]->scale(100.);
      	hr_DstarDplus[i]->scale(100.);
      	hr_DsDstar[i]->scale(100.);
      }

    }

    /// @}


  private:

    /// @name Histograms
    /// @{

    Histo1DGroupPtr _h_pdg411_Dplus_pT_y, _hbr_Dplus;
    Histo1DGroupPtr _h_pdg421_Dzero_pT_y, _hbr_Dzero;
    Histo1DGroupPtr _h_pdg431_Dsplus_pT_y, _hbr_Ds;
    Histo1DGroupPtr _h_pdg413_Dstarplus_pT_y, _hbr_Dstar;
    Estimate1DPtr hr_DplusDzero[5];
    Estimate1DPtr hr_DsDzero[5];
    Estimate1DPtr hr_DstarDzero[5];
    Estimate1DPtr hr_DsDplus[5];
    Estimate1DPtr hr_DstarDplus[5];
    Estimate1DPtr hr_DsDstar[5];

    /// @}

  };


  RIVET_DECLARE_PLUGIN(LHCB_2016_I1490663);

}
