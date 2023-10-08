
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {


  /// @brief Rho meson production in pp and Pb-Pb at 2.76 TeV
  class ALICE_2019_I1672860 : public Analysis {

  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2019_I1672860);

    void init() {

      // Find out the beam type
      const ParticlePair& beam = beams();
      if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON) isHI = false;
      else if (beam.first.pid() == PID::LEAD && beam.second.pid() == PID::LEAD) {
      	isHI = true;
      }
      else {
      	MSG_WARNING("No beam found. You are likely in REENTRANT status.");
      	isHI = true;
      }

      if (isHI) {
        declare(HepMCHeavyIon(), "HepMC");

        declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

        _centrality_regions.clear();
        _centrality_regions = {{0., 20.},   {20., 40.},  {40., 60.}, {60., 80.}};
      }



      // Charged, primary particles with |eta| < 0.5
      declare(ALICE::PrimaryParticles(Cuts::absrap < 0.5 && Cuts::abscharge > 0), "APRIM");

      // Resonances
      declare(UnstableParticles(Cuts::absrap<0.5), "RSN");

      // Booking histograms
      // rho pt spectrum in pp (Table 6 of HEPData and Fig. 5 in Article)
      book(_hist_rho_PP, 6, 1, 1);
      book(_counterSOW_PP, "/TMP/counterSOW_PP");

      // NOTE. pion pt spectrum is needed for the ratio and
      // therefore it requires the same binning as the rho pt spectrum
      std::string name_pion_PP_Fig5 = mkAxisCode(6, 1, 1) + "-pion-pp";
      book(_hist_pion_PP, name_pion_PP_Fig5, refData(6, 1, 1));

      // ratio to pion yield (Table 13, Fig. 9)
      book(_rho_pion_ratio_PP, 13, 1, 1);


      book(_counter_temp, "/TMP/counter.temp"); // counter for PbPb
      book(_counterNcoll_temp, "/TMP/counter.ncoll.temp"); // Ncoll counter for PbPb



      //----------------------------------------------------------------------------------
      // Loop over all histograms
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {

        const string nameCounterPbPb = "/TMP/counter.pbpb." + std::to_string(ihist);
        book(_counterSOW[PBPB][ihist], nameCounterPbPb); // Sum of weights counter for PbPb

        const string nameCounterNcoll = "/TMP/counter.ncoll." + std::to_string(ihist);
        book(_counterNcoll[ihist], nameCounterNcoll); // Ncoll counter for PbPb

        // rho pt spectra in PbPb (Tables 7-10 in HEPData and Fig. 6 in Article)
        book(_hist_rho[PBPB][ihist], ihist+7, 1, 1);

        // NOTE. pion pt spectra are needed for the ratios and
        // therefore require the same binning as the rho pt spectra
        const string name_pion_PbPb = "/TMP/"+mkAxisCode(ihist+7,1,1) + "-pion";
        book(_hist_pion[PBPB][ihist], name_pion_PbPb, refData(ihist+7, 1, 1));

        // NOTE. Only two out of four ratios (0-20% and 60-80%) are presented in the article
        // HEPData (Tables 14, 15) and Article (Fig. 10)
        const string name_rho_pion_ratio_PbPb = "/TMP/"+mkAxisCode(ihist+7,1,1) + "-rho_pion_ratio";
        if (ihist==0)       book(_rho_pion_ratio[PBPB][ihist], 14, 1, 1);
        else if (ihist==3)  book(_rho_pion_ratio[PBPB][ihist], 15, 1, 1);
        else                book(_rho_pion_ratio[PBPB][ihist], name_rho_pion_ratio_PbPb, refData(ihist+7, 1, 1).xEdges());

        // Next lines are for RAA
        // Initialize pp objects. In principle, only one pp histogram would be
        // needed since centrality does not make any difference here. However,
        // in some cases in this analysis the binning differ from each other,
        // so this is easy-to-implement way to account for that.
        const string nameCounterpp = "/TMP/counter.pp." + std::to_string(ihist);
        book(_counterSOW[PP][ihist], nameCounterpp); // Sum of weights counter for pp
        const string namePP = "/TMP/"+mkAxisCode(ihist+7,1,1) + "-pp";
        book(_hist_rho[PP][ihist], namePP, refData(ihist+7, 1, 1));
        const string name_pion_PP = "/TMP/"+mkAxisCode(ihist+7,1,1) + "-pion-pp";
        book(_hist_pion[PP][ihist], name_pion_PP, refData(ihist+7, 1, 1));
        const string name_rho_pion_ratio_PP = "/TMP/"+mkAxisCode(ihist+7,1,1) + "-rho_pion_ratio-pp";
        book(_rho_pion_ratio[PP][ihist], name_rho_pion_ratio_PP, refData(ihist+7, 1, 1).xEdges());
        // RAA are in Tables 16-19 in HEPData and Fig. 11 in Article
        book(_hist_RAA[ihist], ihist+16, 1, 1);
      } // end loop over histograms

      // integrated yields, ratios and mean pt vs. multiplicity
      // Table 11, Fig. 8, left
      book(_hist_integrated_rho_pion_ratio, 11, 1, 1);
      book(_hist_integrated_yield_rho, "/TMP/integrated_yield_rho", refData( 11, 1, 1));
      book(_hist_integrated_yield_pion, "/TMP/integrated_yield_pion", refData( 11, 1, 1));
      book(_hist_mean_pt_rho, 12, 1, 1);

    } // end init

    // Perform the per-event analysis
    void analyze(const Event& event) {

      if (int_edges.empty()) {
        int_edges.push_back(_hist_integrated_yield_pion->bin(1).xMid());
        int_edges.push_back(_hist_integrated_yield_pion->bin(3).xMid());
        int_edges.push_back(_hist_integrated_yield_pion->bin(5).xMid());
        int_edges.push_back(_hist_integrated_yield_pion->bin(7).xMid());
        int_edges.push_back(_hist_integrated_yield_pion->bin(9).xMid());
      }

      // Charged, primary particles in eta range of |eta| < 0.5
      Particles chargedParticles = apply<ALICE::PrimaryParticles>(event,"APRIM").particlesByPt();

      // Resonances
      const UnstableParticles &rsn = apply<UnstableParticles>(event, "RSN");

      if (isHI) {

        const HepMCHeavyIon &hi = apply<HepMCHeavyIon>(event, "HepMC");
        if (!hi.ok()) {
          MSG_WARNING("HEPMC Heavy ion container needed for this analysis, but not "
                "found for this event. Skipping.");
          vetoEvent;
        }

        _counter_temp->fill(hi.Ncoll());
        _counterNcoll_temp->fill(hi.Ncoll());


        // Prepare centrality projection and value
        const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
        double centr = centrProj();
        // Veto event for too large centralities since those are not used
        // in the analysis at all
        if ((centr < 0.) || (centr > 80.)) vetoEvent;

        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {

          const double low_edge_rho = _hist_rho[PBPB][ihist]->xMin();
          const double high_edge_rho = _hist_rho[PBPB][ihist]->xMax();
          const double low_edge_pion = _hist_pion[PBPB][ihist]->xMin();
          const double high_edge_pion = _hist_pion[PBPB][ihist]->xMax();


          if (inRange(centr, _centrality_regions[ihist].first, _centrality_regions[ihist].second)) {

            _counterSOW[PBPB][ihist]->fill();
            _counterNcoll[ihist]->fill(hi.Ncoll());

            for (const Particle &p : rsn.particles()) {
              if (p.abspid() != PID::RHO0)  continue;
              _hist_integrated_yield_rho->fill(int_edges[4-ihist]);
              const double pT = p.pT()/GeV;
              _hist_mean_pt_rho->fill(int_edges[4-ihist], pT);
              if (pT > low_edge_rho && pT < high_edge_rho) {
                _hist_rho[PBPB][ihist]->fill(pT);
              } // condition on pT
            } // end loop over resonances

            //----------------------------------------------------------------------------------
            for (const Particle& p : chargedParticles) {
              if (p.abspid() != PID::PIPLUS)  continue;
              _hist_integrated_yield_pion->fill(int_edges[4-ihist]);
              const double pT = p.pT()/GeV;
              if (pT > low_edge_pion && pT < high_edge_pion) {
                _hist_pion[PBPB][ihist]->fill(pT);
              }
            } // end loop over charged primary particles
          } // centrality
        } // histo loop
      } // end PbPb event
      else { // PP event

        _counterSOW_PP->fill();

        double low_edge_rho = _hist_rho_PP->xMin();
        double high_edge_rho = _hist_rho_PP->xMax();

        for (const Particle &p : rsn.particles()) {
          if (p.abspid() != PID::RHO0)  continue;
          _hist_integrated_yield_rho->fill(int_edges[0]);  // fill first bin for pp
          const double pT = p.pT()/GeV;
          _hist_mean_pt_rho->fill(int_edges[0], pT);
          if (pT > low_edge_rho && pT < high_edge_rho) {
            _hist_rho_PP->fill(pT);
          } // condition on pT
        } // end loop over resonances

        double low_edge_pion = _hist_pion_PP->xMin();
        double high_edge_pion = _hist_pion_PP->xMax();

        for (const Particle& p : chargedParticles) {
          if (p.abspid() != PID::PIPLUS)  continue;
          _hist_integrated_yield_pion->fill(int_edges[0]); // fill first bin for pp
          const double pT = p.pT()/GeV;
          if (pT > low_edge_pion && pT < high_edge_pion) {
            _hist_pion_PP->fill(pT);
          }
        } // end loop over charged primary particles

        // next histograms are needed only for the RAA
        // NOTE. for completeness and consistency checks also histo for pions are provided
        // although they are not needed (see above for rho/pion ratio in pp).
        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {

          _counterSOW[PP][ihist]->fill();

          low_edge_rho = _hist_rho[PP][ihist]->xMin();
          high_edge_rho = _hist_rho[PP][ihist]->xMax();

          for (const Particle &p : rsn.particles()) {
            if (p.abspid() != PID::RHO0)  continue;
            double pT = p.pT()/GeV;
            if (pT > low_edge_rho && pT < high_edge_rho) {
              _hist_rho[PP][ihist]->fill(pT);
            } // condition on pT
          } // end loop over resonances

          low_edge_pion = _hist_pion[PP][ihist]->xMin();
          high_edge_pion = _hist_pion[PP][ihist]->xMax();

          for (const Particle& p : chargedParticles) {
            if (p.abspid() != PID::PIPLUS)  continue;
            const double pT = p.pT()/GeV;
            if (pT > low_edge_pion && pT < high_edge_pion) {
                _hist_pion[PP][ihist]->fill(pT);
            }
          } // end loop over charged primary particles
        } // loop over histos
      } // end pp event
    } // end analyze




    // Normalise histograms etc., after the run
    void finalize() {

      if (_counterSOW_PP->sumW() > 0.) {
        scale(_hist_rho_PP,  1. / _counterSOW_PP->sumW());
        scale(_hist_pion_PP, 1. / _counterSOW_PP->sumW());

        if (_hist_rho_PP->numEntries() > 0 && _hist_pion_PP->numEntries() > 0) {
          divide(_hist_rho_PP, _hist_pion_PP, _rho_pion_ratio_PP);
          scale(_rho_pion_ratio_PP, 2.);
        }

      }

      // Scaling of the histograms with their individual weights.
      for (size_t itype = 0; itype < EVENT_TYPES; ++itype ) {
        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {

          if (_counterSOW[itype][ihist]->sumW() > 0.) {
            scale(_hist_rho[itype][ihist],  1./ _counterSOW[itype][ihist]->sumW());
            scale(_hist_pion[itype][ihist], 1./ _counterSOW[itype][ihist]->sumW());

            if (_hist_rho[itype][ihist]->numEntries() > 0 && _hist_pion[itype][ihist]->numEntries() > 0) {
              divide(_hist_rho[itype][ihist], _hist_pion[itype][ihist], _rho_pion_ratio[itype][ihist]);
              scale(_rho_pion_ratio[itype][ihist], 2.);
            }
          }

        } // end histo loop
      } // PP and PBPB


      // Postprocessing for RAA
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
        // If there are entries in histograms for both beam types
        if (_hist_rho[PP][ihist]->numEntries() > 0 && _hist_rho[PBPB][ihist]->numEntries() > 0) {
          // Initialize and fill R_AA histograms
          divide(_hist_rho[PBPB][ihist], _hist_rho[PP][ihist], _hist_RAA[ihist]);

          // Scale by Ncoll. Unfortunately some generators does not provide
          // Ncoll value (eg. JEWEL), so the following scaling will be done
          // only if there are entries in the counters
          double ncoll = _counterNcoll[ihist]->sumW();

          //PHYSICAL REVIEW C 88, 044909 (2013)
          if(ihist==0) ncoll=1210.85;
          else if(ihist==1) ncoll=438.;
          else if(ihist==2) ncoll=127.7;
          else ncoll=26.7;

          double sow = _counterSOW[PBPB][ihist]->sumW();
          if (ncoll > 1e-6 && sow > 1e-6) {
            scale(_hist_RAA[ihist], 1. / ncoll);
          }

        }
      } // loop over histos


      // Postprocessing for integrated yield vs. multiplicity
      if (_hist_integrated_yield_rho->numEntries() > 0. && _hist_integrated_yield_pion->numEntries() > 0.) {
        divide( _hist_integrated_yield_rho, _hist_integrated_yield_pion, _hist_integrated_rho_pion_ratio);
      }



    } // end finalize
    //=================================================================================


  private:

    bool isHI;
    static const int NHISTOS = 4;
    static const int EVENT_TYPES = 2;
    static const int PP = 0;
    static const int PBPB = 1;


    // pt spectrum for rho in pp@2.76 TeV
    // Table 6 in HEPData and Fig. 5 in Article
    // NOTE. The histo for pions is need for the rho/pion ratio
    // Table 9 in HEPData and Fig. 9 in Article
    CounterPtr _counterSOW_PP;
    Histo1DPtr _hist_rho_PP;
    Histo1DPtr _hist_pion_PP;
    Estimate1DPtr _rho_pion_ratio_PP;

    // pt spectra for rho in PbPb@2.76 TeV in 0-20%, 20-40%, 40-60%, 60-80%
    // Tables 7-10 in HEPData and Fig. 6 in Article
    // NOTE1. For EVENT_TYPES=PP, histos are pt spectra for rho in pp@2.76 TeV
    // with the same binning as the pt spectra in PbPb.
    // They are needed for RAA
    // Tables 16-19 in HEPData and Fig. 11 in Article
    // NOTE2. histos for pions are needed for the rho/pion ratio
    // Tables 14 and 15 in HEPData and Fig. 10 in Article
    Histo1DPtr _hist_rho[EVENT_TYPES][NHISTOS];
    Histo1DPtr _hist_pion[EVENT_TYPES][NHISTOS];
    CounterPtr _counterSOW[EVENT_TYPES][NHISTOS];
    CounterPtr _counterNcoll[NHISTOS];

    CounterPtr _counter_temp;
    CounterPtr _counterNcoll_temp;

    Estimate1DPtr _rho_pion_ratio[EVENT_TYPES][NHISTOS];
    Estimate1DPtr _hist_RAA[NHISTOS];

    // integrated yields vs. multiplicity and mean pT
    // Table 11-12, Fig. 8, left and right
    Histo1DPtr _hist_integrated_yield_rho;
    Histo1DPtr _hist_integrated_yield_pion;
    Estimate1DPtr _hist_integrated_rho_pion_ratio;

    Profile1DPtr _hist_mean_pt_rho;


    std::vector<std::pair<double, double>> _centrality_regions;
    vector<double> int_edges;


  };


  RIVET_DECLARE_PLUGIN(ALICE_2019_I1672860);

}
