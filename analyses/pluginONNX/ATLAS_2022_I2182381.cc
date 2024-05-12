// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Tools/RivetONNXrt.hh"

#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief Search for gluino pair-production in multi-b + MET final states
  class ATLAS_2022_I2182381 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2182381);

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections

      const FinalState fs(Cuts::abseta < 4.9);
      const ChargedFinalState cfs(Cuts::abseta < 2.5 && Cuts::pt > 1*GeV);
      declare(cfs, "chargedTracks");

      // leptons
      // n.b. cut too low to account for smearing
      const DirectFinalState bare_leptons((Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      declare(bare_leptons, "bare_leps");
      const DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      LeptonFinder dressed_leps(bare_leptons, photons, 0.1, Cuts::abseta < 2.6
                                  && Cuts::pt > 15*GeV );

      SmearedParticles smearedelectrons(dressed_leps, Cuts::abspid==PID::ELECTRON 
                                          && !Cuts::absetaIn(1.37,1.52) && Cuts::abseta < 2.47,
                                          PARTICLE_EFF_PERFECT, ELECTRON_SMEAR_ATLAS_RUN2);
    	declare(smearedelectrons, "SmearedElec");
	    SmearedParticles recomuons(dressed_leps, Cuts::abspid==PID::MUON && 
                                  Cuts::pt > 20*GeV && Cuts::abseta < 2.5,
                                  MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2);
      declare(recomuons, "SmearedMuon");

      // Small-R jets
      const FinalState fsj(Cuts::abseta < 4.8);
	    FastJets Sj(fsj, JetAlg::ANTIKT, 0.4, JetMuons::NONE);
      SmearedJets SSj(Sj, JET_SMEAR_ATLAS_RUN2, JET_BTAG_MV2C10_77_2016);
    	declare(SSj, "smearedSjet");

      declare(MissingMom(), "Etmiss");

      // Load the NN
      _nn = getONNX(name());
      // n.b. normalisation values taken from simpleAnalysis config file
      // and inserted into the onnx file metadata
      _nn_norm_means = _nn->retrieve<float>("norm_means");
      _nn_norm_devs = _nn->retrieve<float>("norm_devs");
      MSG_DEBUG(*_nn); 

      // Book counters
      for (const string& s : _cc_srs){
        book(_c[s], "_" + s);
      }
      for (const string& s : _cc_vrs){
        book(_c[s], "_" + s);
      }
      for (const string& s : _cc_crs){
        book(_c[s], "_" + s);
      }
      for (const string& s : _nn_srs){
        book(_c[s], "_" + s);
      }
      for (const string& s : _nn_vrs){
        book(_c[s], "_" + s);
      }
      for (const string& s : _nn_crs){
        book(_c[s], "_" + s);
      }

      /// Cutflows
      // n.b. All regions will apply a 3 bjet cut at some point, but this isn't explicit
      // in some of the cutflow tables. I've assumes it's done with the rest of PreSel
      
      // CC analysis
      book(_f_GbbCC, {"CF_Gbb_0l_B", "CF_Gbb_0l_M", "CF_Gbb_0l_C"});
      book(_f_GtbCC, {"CF_Gtb_0l_B", "CF_Gtb_0l_M", "CF_Gtb_0l_C"});
      book(_f_Gtt0lCC, {"CF_Gtt_0l_B", "CF_Gtt_0l_M1", "CF_Gtt_0l_M2", "CF_Gtt_0l_C"});
      book(_f_Gtt1lCC, {"CF_Gtt_1l_B", "CF_Gtt_1l_M1", "CF_Gtt_1l_M2", "CF_Gtt_1l_C"});
      
      book(_f_Gtt0lCC->binAt("CF_Gtt_0l_B"), "CF_Gtt_0l_B", {"Nlepbase=0","dPhi4jmin>0.4","NJet>=5","ETMiss>=600","Meff>=2900","mbjets_Tmin>120","TotJetMass>=300"});
      book(_f_Gtt0lCC->binAt("CF_Gtt_0l_M1"), "CF_Gtt_0l_M1", {"Nlepbase=0","dPhi4jmin>0.4","NJet>=9&&Nbs>=3","ETMiss>=600","Meff>=1700","mbjets_Tmin>120","TotJetMass>=300"});
      book(_f_Gtt0lCC->binAt("CF_Gtt_0l_M2"), "CF_Gtt_0l_M2", {"Nlepbase=0","dPhi4jmin>0.4","NJet>=10&&Nbs>=3","ETMiss>=500","Meff>=1100","mbjets_Tmin>120","TotJetMass>=200"}); 
      book(_f_Gtt0lCC->binAt("CF_Gtt_0l_C"), "CF_Gtt_0l_C", {"Nlepbase=0","dPhi4jmin>0.4","NJets>=10","NBjets>=4","ETMiss>=400","Meff>=800","mbjets_Tmin>180","TotJetMass>=100"});

      book(_f_GbbCC->binAt("CF_Gbb_0l_B"), "CF_Gbb_0l_B", {"Nlepbase=0","dPhi4jmin>0.4","mbjets_Tmin>=130","ETMiss>=550","pTjet>=65","Meff>=2600"});
      book(_f_GbbCC->binAt("CF_Gbb_0l_M"), "CF_Gbb_0l_M", {"Nlepbase=0","dPhi4jmin>0.4","mbjets_Tmin>=130","ETMiss>=550","Meff>=2000"});
      // n.b. Auxiliary table 4 clearly labels Gbb_C meff and met cuts the wrong way round, should be consistent with Gbb_B, Gbb_M
      book(_f_GbbCC->binAt("CF_Gbb_0l_C"), "CF_Gbb_0l_C", {"Nlepbase=0","dPhi4,jmin>0.4","mbjets_Tmin>=130","ETMiss>=550","Meff>=1600"}); 

      book(_f_GtbCC->binAt("CF_Gtb_0l_B"), "CF_Gtb_0l_B", {"ETmissTrigger", "Nlepbase=0", "dPhi4jmin>0.4", "mbjets_Tmin>=130", "Meff>=2500", "ETMiss>=550", "totJetMass>=200"});
      book(_f_GtbCC->binAt("CF_Gtb_0l_M"), "CF_Gtb_0l_M", {"ETmissTrigger", "Nlepbase=0", "dPhi4jmin>0.4", "mbjets_Tmin>=130", "nJets>=6", "nbs>=4", "Meff>=2000", "ETMiss>=550","totJetMass>=200"});
      book(_f_GtbCC->binAt("CF_Gtb_0l_C"), "CF_Gtb_0l_C", {"ETmissTrigger", "Nlepbase=0", "dPhi4jmin>0.4", "mbjets_Tmin>=130", "nJets>=7", "nbs>=4", "Meff>=1300", "ETMiss>=500","totJetMass>=50"}); 

      book(_f_Gtt1lCC->binAt("CF_Gtt_1l_B"), "CF_Gtt_1l_B", {"Nsiglep=1","NJet>=4&&nbs>=3","ETMiss>=600","Meff>=2300","m_trans>=150","mbjets_Tmin>120","TotJetMass>=200"});
      book(_f_Gtt1lCC->binAt("CF_Gtt_1l_M1"), "CF_Gtt_1l_M1", {"Nsiglep=1","NJet>=5&&nbs>=3","ETMiss>=600","Meff>=2000","m_trans>=200","mbjets_Tmin>120","TotJetMass>=200"});
      book(_f_Gtt1lCC->binAt("CF_Gtt_1l_M2"), "CF_Gtt_1l_M2", {"Nsiglep=1","NJet>=8&&nbs>=3","ETMiss>=500","Meff>=1100","m_trans>=200","mbjets_Tmin>120","TotJetMass>=100"}); 
      book(_f_Gtt1lCC->binAt("CF_Gtt_1l_C"), "CF_Gtt_1l_C", {"Nsiglep=1","NJet>=9&&nbs>=3","ETMiss>=300","Meff>=800","m_trans>=150","mbjets_Tmin>120"});

      // NN regions
      book(_f_GbbNN, {"CF_Gbb_2800_1400", "CF_Gbb_2300_1000", "CF_Gbb_2100_1600", "CF_Gbb_2000_1800"});
      book(_f_GttNN, {"CF_Gtt_2100_1", "CF_Gtt_1800_1", "CF_Gtt_2300_1200", "CF_Gtt_1900_1400"});
      const static vector<string> NN_Gbb_cuts{"NBaseLep=0", "dPhi4jmin", "P(Gbb)"};
      const static vector<string> NN_Gtt_cuts{"Common", "P(Gtt)"};
      book(_f_GttNN->binAt("CF_Gtt_2100_1"), "CF_Gtt_2100_1", NN_Gtt_cuts);
      book(_f_GttNN->binAt("CF_Gtt_1800_1"), "CF_Gtt_1800_1", NN_Gtt_cuts);
      book(_f_GttNN->binAt("CF_Gtt_2300_1200"), "CF_Gtt_2300_1200", NN_Gtt_cuts);
      book(_f_GttNN->binAt("CF_Gtt_1900_1400"), "CF_Gtt_1900_1400", NN_Gtt_cuts);
      book(_f_GbbNN->binAt("CF_Gbb_2800_1400"), "CF_Gbb_2800_1400", NN_Gbb_cuts);
      book(_f_GbbNN->binAt("CF_Gbb_2300_1000"), "CF_Gbb_2300_1000", NN_Gbb_cuts);
      book(_f_GbbNN->binAt("CF_Gbb_2100_1600"), "CF_Gbb_2100_1600", NN_Gbb_cuts);
      book(_f_GbbNN->binAt("CF_Gbb_2000_1800"), "CF_Gbb_2000_1800", NN_Gbb_cuts);
     
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get (smallR) jets
      Jets smallJets = apply<JetFinder>(event, "smearedSjet").jetsByPt(Cuts::pt > 30*GeV && Cuts::abseta < 2.8);

      // Let's have a go at applying a JVT like-cut.
      idiscard(smallJets, [](const Jet &j){return rand01() > ATLAS_JVT_EFF_MEDIUM(j);});

      // Get the Leptons.
      Particles smeared_muons = apply<ParticleFinder>(event, "SmearedMuon").particles();
      // Manually apply electron efficiencies so I can get seperate but consistent tight/loose definitions
      // TODO: Can I think of a projection to do this?
      Particles true_electrons = apply<ParticleFinder>(event, "SmearedElec").particles(Cuts::abspid == PID::ELECTRON
                                                                            && Cuts::abseta < 2.47 
                                                                            && !Cuts::absetaIn(1.37, 1.52)
                                                                            && Cuts::pt > 20*GeV);
      Particles tight_electrons, loose_electrons;
      {
        for (const Particle & e : true_electrons){
          const double rnd = rand01();
          if (ELECTRON_EFF_ATLAS_RUN2_TIGHT(e) > rnd){
            tight_electrons.push_back(e);
          }
          if (ELECTRON_EFF_ATLAS_RUN2_LOOSE(e) > rnd){
            loose_electrons.push_back(e);
          }
        }
      }

      // Overlap removal:
      // jets within 0.2 of an electron (unless btagged)
      // TODO: technically I should be using a "looser" definition of btagged for this...??
      idiscardIfAny(smallJets, loose_electrons, [](const Jet & j, const Particle &e)
        {return deltaR(j,e) < 0.2 && !j.bTagged();});
      // remaining electrons within 0.4 of a jet
      idiscardIfAnyDeltaRLess(loose_electrons, smallJets, 0.4);
      idiscardIfAnyDeltaRLess(tight_electrons, smallJets, 0.4);

      // muons within 0.4 of a jet, if the jet has > 2 consituents
      idiscardIfAny(smeared_muons, smallJets, [](const Particle &mu, const Jet & j)
        {return deltaR(j,mu) < 0.4 && j.size() > 2;});
      // jets within 0.4 of a muon if they have 2 or fewer constituents
      idiscardIfAny(smallJets, smeared_muons, [](const Jet & j, const Particle &mu)
        {return deltaR(j,mu) < 0.4 && j.size() <= 2;});

      // bjets
      const Jets bjets = sortByPt(select(select(smallJets, Cuts::abseta < 2.5), hasBTag()));
      // Signal lepton track based requirements - Charged tracks final state required:
      Particles cfs = apply<ParticleFinder>(event, "chargedTracks").particles(Cuts::pt > 1*GeV && Cuts::abseta < 2.5);
      // Signal Muons: FixedCutTightTrackOnly - definition inferred from arXiv:1603.05598 table 2
      Particles tight_muons = discard(smeared_muons, [&cfs](const Particle& mu){
        const double muRadius = min(0.3, (10*GeV)/mu.pt());
        return (1.06*mu.pt() < std::accumulate(cfs.begin(), cfs.end(), 0.0, 
          [&mu, &muRadius](const double tot, const Particle & track){
            return deltaR(mu, track) < muRadius ? tot+track.pt() : tot ;
        }));
      });
      // Signal electrons: Fix(loose) - from 1902.04655 table 4.
      // n.b. only pT (not Et) part of the criteria implemented - How to define the Et cone in rivet I have no idea
      // Magic 0.2 given in text of paper.
      idiscard(tight_electrons, [&cfs](const Particle& el){
        const double elRadius = min(0.2, (10*GeV)/el.pt());
        return (1.15*el.pt() < std::accumulate(cfs.begin(), cfs.end(), 0.0, 
          [&el, &elRadius](const double tot, const Particle & track){
            return (deltaR(el, track) < elRadius && (deltaPhi(el, track) > 0.1 && deltaEta(el, track) > 0.05)) ? tot+track.pt() : tot ;
        }));
      });

      //Recluster to make large-R jets
      fastjet::JetDefinition jdef(fastjet::antikt_algorithm, 0.8);
      ClusterSequence cseq(smallJets, jdef);
      PseudoJets largeRPseudoJets = cseq.inclusive_jets();
      Jets largeRjets = trimJetsFrac(largeRPseudoJets, 0.1);
      iselect(largeRjets, Cuts::abseta < 2.0 && Cuts::pt > 100*GeV);
      isortByPt(largeRjets);
      const size_t nlargeRjets = largeRjets.size();

      const MissingMomentum ETMiss = apply<MissingMomentum>("Etmiss", event);
      const double etmiss = ETMiss.scalarPtMiss();

      const size_t nSigLeptons = tight_muons.size()+tight_electrons.size();
      const size_t nBaseLeptons = smeared_muons.size()+loose_electrons.size();
      const Particles sigleptons = sortByPt(tight_muons + tight_electrons);
      const Particles baseleptons = sortByPt(smeared_muons + loose_electrons);

      bool zjcr = false; // Does the event qualify for a zjcr
      bool onlyzjcr = false; // Does the event qualify for zjcr but would otherwise fail preselection
      ThreeMomentum ETMisshat{};
      // A small (and hopefully quick) aside: Gbb NN analysis has four very weird CRs:
      if (nSigLeptons == 2 &&
          sigleptons[0].pid() == -sigleptons[1].pid()){
        const double combined_lepton_mass = (sigleptons[0].mom()+sigleptons[1].mom()).mass();
        const double combined_lepton_pt = (sigleptons[0].mom()+sigleptons[1].mom()).pt();
        if(combined_lepton_mass > 60*GeV && combined_lepton_mass < 120*GeV
        && combined_lepton_pt > 70*GeV ){
          ETMisshat = ETMiss.vectorEtMiss() + sigleptons[0].p3() + sigleptons[1].p3();
          if (ETMisshat.pt() > 200*GeV ){
            zjcr = true;
            if (etmiss < 200){
              onlyzjcr = true;
            }
          }
        }
      }
      

      // Universal preselection
      // At least four jets
      // (TODO except for Zjets CR?)
      if (smallJets.size() < 4) vetoEvent;
      // At least three of which btagged:
      if (bjets.size() < 3) vetoEvent;
      // ETmiss > 200GeV
      if ((etmiss <= 200*GeV) && !onlyzjcr) vetoEvent;

      if (!onlyzjcr){
        _f_GbbNN->groupfillinit(); _f_GttNN->groupfillinit();
        _f_GbbCC->groupfillinit(); _f_Gtt0lCC->groupfillinit(); _f_Gtt1lCC->groupfillinit(); _f_GtbCC->groupfillinit();

        // Gtb 0l 1st cutflow entry
        _f_GtbCC->groupfillnext();
      }

      // Define some important variables
      const double meff = (etmiss +
        accumulate(smallJets.begin(), smallJets.end(), 0., [](double tot, const Jet &j){return tot+j.pt();}) +
        accumulate(tight_electrons.begin(), tight_electrons.end(), 0., [](double tot, const Particle &e){return tot+e.pt();}) +
        accumulate(smeared_muons.begin(), smeared_muons.end(), 0.,[](double tot, const Particle &mu){return tot+mu.pt();}));

      const double totJetMass = std::accumulate(largeRjets.begin(), largeRjets.end(), 0.0, 
                                  [](double tot, const Jet &j){return tot + j.mass();});

      const size_t nbs = bjets.size();

      const size_t njets = smallJets.size();

      const double mBjets_Tmin = min({mT(ETMiss.vectorEtMiss(), bjets[0].p3()),
                                      mT(ETMiss.vectorEtMiss(), bjets[1].p3()),
                                      mT(ETMiss.vectorEtMiss(), bjets[2].p3())});
      
      /////////////////////////////////////////////////////////////////////////
      // NN Analysis
      /////////////////////////////////////////////////////////////////////////
      // The NN analysis has a more inclusive selection -> let's do it first
      {
        // We need to pass criteria for either the 0l OR 1l regions
        const double dPhi4jmin = min({deltaPhi(ETMiss.vectorMET(), smallJets[0]),
                                  deltaPhi(ETMiss.vectorMET(), smallJets[1]),
                                  deltaPhi(ETMiss.vectorMET(), smallJets[2]),
                                  deltaPhi(ETMiss.vectorMET(), smallJets[3])});


        // TODO: what is the right balance of width vs height vs organised rationally for this massive vector?
        vector<float> nn_input{
          // small jet variables.
          (float)smallJets[0].pt(), (float)smallJets[0].eta(), (float)smallJets[0].phi(MINUSPI_PLUSPI), (float)smallJets[0].mass(),  (float)smallJets[0].bTagged(),
          (float)smallJets[1].pt(), (float)smallJets[1].eta(), (float)smallJets[1].phi(MINUSPI_PLUSPI), (float)smallJets[1].mass(),  (float)smallJets[1].bTagged(),
          (float)smallJets[2].pt(), (float)smallJets[2].eta(), (float)smallJets[2].phi(MINUSPI_PLUSPI), (float)smallJets[2].mass(),  (float)smallJets[2].bTagged(),
          (float)smallJets[3].pt(), (float)smallJets[3].eta(), (float)smallJets[3].phi(MINUSPI_PLUSPI), (float)smallJets[3].mass(),  (float)smallJets[3].bTagged(),
          (njets > 4 ? (float)smallJets[4].pt() : 0.0f), njets > 4 ? (float)smallJets[4].eta() : 0.f, njets > 4 ? (float)smallJets[4].phi(MINUSPI_PLUSPI) : 0.f, njets > 4 ? (float)smallJets[4].mass() : 0.f,  njets > 4 ? (float)smallJets[4].bTagged() : 0.f,
          (njets > 5 ? (float)smallJets[5].pt() : 0.0f), njets > 5 ? (float)smallJets[5].eta() : 0.f, njets > 5 ? (float)smallJets[5].phi(MINUSPI_PLUSPI) : 0.f, njets > 5 ? (float)smallJets[5].mass() : 0.f,  njets > 5 ? (float)smallJets[5].bTagged() : 0.f,
          (njets > 6 ? (float)smallJets[6].pt() : 0.0f), njets > 6 ? (float)smallJets[6].eta() : 0.f, njets > 6 ? (float)smallJets[6].phi(MINUSPI_PLUSPI) : 0.f, njets > 6 ? (float)smallJets[6].mass() : 0.f,  njets > 6 ? (float)smallJets[6].bTagged() : 0.f,
          (njets > 7 ? (float)smallJets[7].pt() : 0.0f), njets > 7 ? (float)smallJets[7].eta() : 0.f, njets > 7 ? (float)smallJets[7].phi(MINUSPI_PLUSPI) : 0.f, njets > 7 ? (float)smallJets[7].mass() : 0.f,  njets > 7 ? (float)smallJets[7].bTagged() : 0.f,
          (njets > 8 ? (float)smallJets[8].pt() : 0.0f), njets > 8 ? (float)smallJets[8].eta() : 0.f, njets > 8 ? (float)smallJets[8].phi(MINUSPI_PLUSPI) : 0.f, njets > 8 ? (float)smallJets[8].mass() : 0.f,  njets > 8 ? (float)smallJets[8].bTagged() : 0.f,
          (njets > 9 ? (float)smallJets[9].pt() : 0.0f), njets > 9 ? (float)smallJets[9].eta() : 0.f, njets > 9 ? (float)smallJets[9].phi(MINUSPI_PLUSPI) : 0.f, njets > 9 ? (float)smallJets[9].mass() : 0.f,  njets > 9 ? (float)smallJets[9].bTagged() : 0.f,
          
          // large Jet variables
          nlargeRjets > 0 ? (float)largeRjets[0].pt() : 0.f, nlargeRjets > 0 ? (float)largeRjets[0].eta() : 0.f, nlargeRjets > 0 ? (float)largeRjets[0].phi(MINUSPI_PLUSPI) : 0.f, nlargeRjets > 0 ? (float)largeRjets[0].mass() : 0.f,
          nlargeRjets > 1 ? (float)largeRjets[1].pt() : 0.f, nlargeRjets > 1 ? (float)largeRjets[1].eta() : 0.f, nlargeRjets > 1 ? (float)largeRjets[1].phi(MINUSPI_PLUSPI) : 0.f, nlargeRjets > 1 ? (float)largeRjets[1].mass() : 0.f,
          nlargeRjets > 2 ? (float)largeRjets[2].pt() : 0.f, nlargeRjets > 2 ? (float)largeRjets[2].eta() : 0.f, nlargeRjets > 2 ? (float)largeRjets[2].phi(MINUSPI_PLUSPI) : 0.f, nlargeRjets > 2 ? (float)largeRjets[2].mass() : 0.f,
          nlargeRjets > 3 ? (float)largeRjets[3].pt() : 0.f, nlargeRjets > 3 ? (float)largeRjets[3].eta() : 0.f, nlargeRjets > 3 ? (float)largeRjets[3].phi(MINUSPI_PLUSPI) : 0.f, nlargeRjets > 3 ? (float)largeRjets[3].mass() : 0.f,

          // lepton variables
          // n.b. using signal leptons following the simpleAnalysis
          nSigLeptons > 0 ? (float)sigleptons[0].pt() : 0.f, nSigLeptons > 0 ? (float)sigleptons[0].eta() : 0.f, nSigLeptons > 0 ? (float)sigleptons[0].phi(MINUSPI_PLUSPI) : 0.f, nSigLeptons > 0 ? (float)sigleptons[0].mass() : 0.f,
          nSigLeptons > 1 ? (float)sigleptons[1].pt() : 0.f, nSigLeptons > 1 ? (float)sigleptons[1].eta() : 0.f, nSigLeptons > 1 ? (float)sigleptons[1].phi(MINUSPI_PLUSPI) : 0.f, nSigLeptons > 1 ? (float)sigleptons[1].mass() : 0.f,
          nSigLeptons > 2 ? (float)sigleptons[2].pt() : 0.f, nSigLeptons > 2 ? (float)sigleptons[2].eta() : 0.f, nSigLeptons > 2 ? (float)sigleptons[2].phi(MINUSPI_PLUSPI) : 0.f, nSigLeptons > 2 ? (float)sigleptons[2].mass() : 0.f,
          nSigLeptons > 3 ? (float)sigleptons[3].pt() : 0.f, nSigLeptons > 3 ? (float)sigleptons[3].eta() : 0.f, nSigLeptons > 3 ? (float)sigleptons[3].phi(MINUSPI_PLUSPI) : 0.f, nSigLeptons > 3 ? (float)sigleptons[3].mass() : 0.f,

          // MET variables
          (float)etmiss, (float)ETMiss.vectorEtMiss().phi(MINUSPI_PLUSPI),

          // Params:
          0.,0.,0.
        };

        // Get outputs 
        normalise_kinematic_nn_input(nn_input);
        map<string, vector<float>> nn_outputs;
        for (const pair<size_t, size_t> & params : _gtt_nn_params){
          nn_input[84] = 1.0; nn_input[85] = params.first; nn_input[86] = params.second;
          normalise_parameter_nn_input(nn_input);
          nn_outputs["Gtt_"+to_string(params.first)+"_"+to_string(params.second)] = _nn->compute(nn_input);
        }
        for (const pair<size_t, size_t> & params : _gbb_nn_params){
          nn_input[84] = 0.0; nn_input[85] = params.first; nn_input[86] = params.second;
          normalise_parameter_nn_input(nn_input);
          nn_outputs["Gbb_"+to_string(params.first)+"_"+to_string(params.second)] = _nn->compute(nn_input);
        }
        // Repeat again for zjcr with leptons zeroed and using etmisshat
        if (zjcr){
          for (const pair<size_t, size_t> & params : _gbb_nn_params){
            // Parameters first
            nn_input[84] = 0.0; nn_input[85] = params.first; nn_input[86] = params.second;
            normalise_parameter_nn_input(nn_input);
            // Zero leptons, set etmiss to etmisshat and then renorm
            zero_lepton_kinematic_nn_input(nn_input);
            nn_input[82] = ETMisshat.pt();
            nn_input[83] = ETMisshat.phi(MINUSPI_PLUSPI);
            normalise_lepmet_kinematic_nn_input(nn_input);
            nn_outputs["Gbb_"+to_string(params.first)+"_"+to_string(params.second)+"_CRZscore"] = _nn->compute(nn_input);
          }
        }

        // Fill Gtt NN sig, control and validation regions
        if (!onlyzjcr){
          if (nSigLeptons >= 1 || (dPhi4jmin >= 0.4 && nBaseLeptons == 0)){
            _f_GttNN->groupfillnext();
            if(_f_GttNN->fillnext("CF_Gtt_2100_1", nn_outputs["Gtt_2100_1"][0] > 0.9997)){
              _c["SR_Gtt_2100_1"]->fill();
            }
            else if (nn_outputs["Gtt_2100_1"][0] > 0.68 && nn_outputs["Gtt_2100_1"][0] < 0.86 &&
                      log10(nn_outputs["Gtt_2100_1"][0]) >= -1.8 && meff >= 2000 ){
              _c["CR_Gtt_2100_1"]->fill();
            }
            else if (nn_outputs["Gtt_2100_1"][0] > 0.86 && nn_outputs["Gtt_2100_1"][0] < 0.9997
                      && meff >= 2000 ){
              _c["VR_Gtt_2100_1"]->fill();
            }

            if(_f_GttNN->fillnext("CF_Gtt_1800_1", nn_outputs["Gtt_1800_1"][0] > 0.9997)){
              _c["SR_Gtt_1800_1"]->fill();
            }
            else if (nn_outputs["Gtt_1800_1"][0] > 0.73 && nn_outputs["Gtt_1800_1"][0] < 0.89 &&
                      log10(nn_outputs["Gtt_1800_1"][0]) >= -2.0 && meff >= 2000 ){
              _c["CR_Gtt_1800_1"]->fill();
            }
            else if (nn_outputs["Gtt_1800_1"][0] > 0.89 && nn_outputs["Gtt_1800_1"][0] < 0.9997
                      && meff >= 2000 ){
              _c["VR_Gtt_1800_1"]->fill();
            }

            if(_f_GttNN->fillnext("CF_Gtt_2300_1200", nn_outputs["Gtt_2300_1200"][0] > 0.9993)){
              _c["SR_Gtt_2300_1200"]->fill();
            }
            else if (nn_outputs["Gtt_2300_1200"][0] > 0.78 && nn_outputs["Gtt_2300_1200"][0] < 0.83 &&
                      log10(nn_outputs["Gtt_2300_1200"][0]) >= -1.6 && meff >= 1400 ){
              _c["CR_Gtt_2300_1200"]->fill();
            }
            else if (nn_outputs["Gtt_2300_1200"][0] > 0.83 && nn_outputs["Gtt_2300_1200"][0] < 0.9993
                      && meff >= 1800 ){
              _c["VR_Gtt_2300_1200"]->fill();
            }

            if(_f_GttNN->fillnext("CF_Gtt_1900_1400", nn_outputs["Gtt_1900_1400"][0] > 0.9987)){
              _c["SR_Gtt_1900_1400"]->fill();
            }
            else if (nn_outputs["Gtt_1900_1400"][0] > 0.78 && nn_outputs["Gtt_1900_1400"][0] < 0.8 &&
                      log10(nn_outputs["Gtt_1900_1400"][0]) >= -1.4 && meff >= 800  && totJetMass < 700){
              _c["CR_Gtt_1900_1400"]->fill();
            }
            else if (nn_outputs["Gtt_1900_1400"][0] > 0.8 && nn_outputs["Gtt_1900_1400"][0] < 0.9987
                      && meff >= 800  && totJetMass < 700){
              _c["VR_Gtt_1900_1400"]->fill();
            }
          }

          // Fill  Gbb NN sig,CR & VR regions
          if (nBaseLeptons ==  0 ){
            _f_GbbNN->groupfillnext();
            if (_f_GbbNN->fillnext("CF_Gbb_2800_1400", {dPhi4jmin >= 0.6, nn_outputs["Gbb_2800_1400"][1] > 0.999})){
              _c["SR_Gbb_2800_1400"]->fill();
            }
            else if (dPhi4jmin >= 0.5 && nn_outputs["Gbb_2800_1400"][1] > 0.43 && nn_outputs["Gbb_2800_1400"][1] < 0.76 
                      && log10(nn_outputs["Gbb_2800_1400"][2]) >= -0.7 && meff >= 1400*GeV && totJetMass < 800*GeV){
              _c["CR_Gbb_2800_1400"]->fill();
            }
            else if (dPhi4jmin >= 0.5 && nn_outputs["Gbb_2800_1400"][1] > 0.76 && nn_outputs["Gbb_2800_1400"][1] < 0.999 
                      && log10(nn_outputs["Gbb_2800_1400"][6]) < -1.7 && meff >= 2500*GeV && totJetMass < 800*GeV){
              _c["VR1_Gbb_2800_1400"]->fill();
            }
            else if (dPhi4jmin >= 0.5 && nn_outputs["Gbb_2800_1400"][1] > 0.76 && nn_outputs["Gbb_2800_1400"][1] < 0.999 
                      && log10(nn_outputs["Gbb_2800_1400"][6]) >= -1.7){
              _c["VR2_Gbb_2800_1400"]->fill();
            }

            if (_f_GbbNN->fillnext("CF_Gbb_2300_1000", {dPhi4jmin >= 0.6, nn_outputs["Gbb_2300_1000"][1] > 0.9994})){
              _c["SR_Gbb_2300_1000"]->fill();
            }
            else if (dPhi4jmin >= 0.5 && nn_outputs["Gbb_2300_1000"][1] > 0.52 && nn_outputs["Gbb_2300_1000"][1] < 0.77 
                      && log10(nn_outputs["Gbb_2300_1000"][2]) >= -0.8 && meff >= 1400*GeV && totJetMass < 800*GeV){
              _c["CR_Gbb_2300_1000"]->fill();
            }
            else if (dPhi4jmin >= 0.5 && nn_outputs["Gbb_2300_1000"][1] > 0.77 && nn_outputs["Gbb_2300_1000"][1] < 0.9994
                      && log10(nn_outputs["Gbb_2300_1000"][6]) < -1.3 && meff >= 2400*GeV && totJetMass < 800*GeV){
              _c["VR1_Gbb_2300_1000"]->fill();
            }
            else if (dPhi4jmin >= 0.5 && nn_outputs["Gbb_2300_1000"][1] > 0.77 && nn_outputs["Gbb_2300_1000"][1] < 0.9994
                      && log10(nn_outputs["Gbb_2300_1000"][6]) >= -1.3){
              _c["VR2_Gbb_2300_1000"]->fill();
            }

            //N.b., because we want to be efficient with our analysis, the dPhi4jmin cut has already been made
            if (_f_GbbNN->fillnext("CF_Gbb_2100_1600", {dPhi4jmin >= 0.4, nn_outputs["Gbb_2100_1600"][1] > 0.9993})){
              _c["SR_Gbb_2100_1600"]->fill();
            }
            else if (dPhi4jmin >= 0.4 && nn_outputs["Gbb_2100_1600"][1] > 0.88 && nn_outputs["Gbb_2100_1600"][1] < 0.91 
                      && log10(nn_outputs["Gbb_2100_1600"][2]) >= -1.3 && meff >= 800*GeV && totJetMass < 500*GeV){
              _c["CR_Gbb_2100_1600"]->fill();
            }
            else if (dPhi4jmin >= 0.4 && nn_outputs["Gbb_2100_1600"][1] > 0.91 && nn_outputs["Gbb_2100_1600"][1] < 0.9993
                      && log10(nn_outputs["Gbb_2100_1600"][6]) < -1.4 && meff >= 800*GeV && totJetMass < 500*GeV){
              _c["VR1_Gbb_2100_1600"]->fill();
            }
            else if (dPhi4jmin >= 0.4 && nn_outputs["Gbb_2100_1600"][1] > 0.91 && nn_outputs["Gbb_2100_1600"][1] < 0.9993
                      && log10(nn_outputs["Gbb_2100_1600"][6]) >= -1.4){
              _c["VR2_Gbb_2100_1600"]->fill();
            }

            //N.b., because we want to be efficient with our analysis, the dPhi4jmin cut has already been made
            if (_f_GbbNN->fillnext("CF_Gbb_2000_1800", {dPhi4jmin >= 0.4, nn_outputs["Gbb_2000_1800"][1] > 0.997})){
              _c["SR_Gbb_2000_1800"]->fill();
            }
            else if (dPhi4jmin >= 0.4 && nn_outputs["Gbb_2000_1800"][1] > 0.92 && nn_outputs["Gbb_2000_1800"][1] < 0.93
                      && log10(nn_outputs["Gbb_2000_1800"][2]) >= -1.9 && meff >= 400*GeV && totJetMass < 400*GeV){
              _c["CR_Gbb_2000_1800"]->fill();
            }
            else if (dPhi4jmin >= 0.4 && nn_outputs["Gbb_2000_1800"][1] > 0.93 && nn_outputs["Gbb_2000_1800"][1] < 0.9997
                      && log10(nn_outputs["Gbb_2000_1800"][6]) < -1.4 && meff >= 400*GeV && totJetMass < 400*GeV){
              _c["VR1_Gbb_2000_1800"]->fill();
            }
            else if (dPhi4jmin >= 0.4 && nn_outputs["Gbb_2000_1800"][1] > 0.93 && nn_outputs["Gbb_2000_1800"][1] < 0.9997
                      && log10(nn_outputs["Gbb_2000_1800"][6]) >= -1.4){
              _c["VR2_Gbb_2000_1800"]->fill();
            }
          }
        }

        // Fill the Gbb NN Z CRs
        if (zjcr){
          if (nn_outputs["Gbb_2000_1800_CRZscore"][1] >= 0.9 &&
               log10(nn_outputs["Gbb_2000_1800_CRZscore"][2]) < -2.2 &&
                log10(nn_outputs["Gbb_2000_1800_CRZscore"][6]) >= -3.7 ){
            _c["CRZ_Gbb_2000_1800"]->fill();
          }
          if (nn_outputs["Gbb_2100_1600_CRZscore"][1] >= 0.9 &&
               log10(nn_outputs["Gbb_2100_1600_CRZscore"][2]) < -1.7 &&
                log10(nn_outputs["Gbb_2100_1600_CRZscore"][6]) >= -4.9 ){
            _c["CRZ_Gbb_2100_1600"]->fill();
          }
          if (nn_outputs["Gbb_2300_1000_CRZscore"][1] >= 0.8 &&
               log10(nn_outputs["Gbb_2300_1000_CRZscore"][2]) < -1.3 &&
                log10(nn_outputs["Gbb_2300_1000_CRZscore"][6]) >= -2.7 ){
            _c["CRZ_Gbb_2300_1000"]->fill();
          }
          if (nn_outputs["Gbb_2800_1400_CRZscore"][1] >= 0.6 &&
               log10(nn_outputs["Gbb_2800_1400_CRZscore"][2]) < -1.3 &&
                log10(nn_outputs["Gbb_2800_1400_CRZscore"][6]) >= -3.0 ){
            _c["CRZ_Gbb_2800_1400"]->fill();
          }
        }
      }
      // The special zjcr events are no longer needed
      if (onlyzjcr) vetoEvent;

      /////////////////////////////////////////////////////////////////////////
      // CC Analysis
      /////////////////////////////////////////////////////////////////////////

      // CC 0l channel
      // n.b. simpleAnalysis uses sigLeps here, but that's not what the paper says.
      // Paper text says base leptons, as do tables 3, 5 and 6, and supplmentary tables 4 & 5.
      // But suplementary table 2 (Gtt cutflows) says it should be siglep.
      if (nBaseLeptons == 0){
        // Fill 0L cutflow
        _f_Gtt0lCC->groupfillnext(); _f_GbbCC->groupfillnext(); _f_GtbCC->groupfillnext();

        // Require dPhi4jmin > 0.4;
        const double dPhi4jmin = min({deltaPhi(ETMiss.vectorMET(), smallJets[0]),
                                  deltaPhi(ETMiss.vectorMET(), smallJets[1]),
                                  deltaPhi(ETMiss.vectorMET(), smallJets[2]),
                                  deltaPhi(ETMiss.vectorMET(), smallJets[3])});
        if (dPhi4jmin  < 0.4) vetoEvent;

        // Fill dPhi4jmin cutflow
        // Gtt
        _f_Gtt0lCC->groupfillnext(); _f_GbbCC->groupfillnext(); _f_GtbCC->groupfillnext();

        // Start filling SRs
        // CC Gtt 0l regions
        if (_f_Gtt0lCC->fillnext("CF_Gtt_0l_B", {njets >= 5, etmiss > 600*GeV, meff >= 2900*GeV, mBjets_Tmin > 120*GeV, totJetMass >= 300*GeV})){
          _c["SR_Gtt_0l_B"]->fill();
        } else if (njets >= 5 && etmiss > 250*GeV && meff > 2000*GeV && totJetMass < 300*GeV){
          _c["VR_Gtt_0l_B"]->fill();
        }
        if (_f_Gtt0lCC->fillnext("CF_Gtt_0l_M1", {njets >= 9 && nbs >= 3, etmiss > 600*GeV, meff >= 1700*GeV, mBjets_Tmin > 120*GeV, totJetMass >= 300*GeV})){
          _c["SR_Gtt_0l_M1"]->fill();
        } else if (njets >= 9 && nbs >= 3 && etmiss > 300*GeV && meff > 1400*GeV && totJetMass < 300*GeV){
          _c["VR_Gtt_0l_M1"]->fill();
        }
        if (_f_Gtt0lCC->fillnext("CF_Gtt_0l_M2", {njets >= 10 && nbs >= 3, etmiss > 500*GeV, meff >= 1100*GeV, mBjets_Tmin > 120*GeV, totJetMass >= 200*GeV})){
          _c["SR_Gtt_0l_M2"]->fill();
        } else if (njets >= 10 && nbs >= 3 && etmiss > 300*GeV && meff >= 800*GeV && totJetMass < 200*GeV){
          _c["VR_Gtt_0l_M2"]->fill();
        }
        if (_f_Gtt0lCC->fillnext("CF_Gtt_0l_C", {njets >= 10, nbs >= 4, etmiss > 400*GeV, meff >= 800*GeV, mBjets_Tmin > 180*GeV, totJetMass >= 100*GeV})){
          _c["SR_Gtt_0l_C"]->fill();
        } else if (njets >= 10 && nbs >= 4 && etmiss > 200*GeV && meff >= 800*GeV && totJetMass < 100*GeV){
          _c["VR_Gtt_0l_M2"]->fill();
        }
        // Gbb 0L
        // CC Gbb 0l SR & VRs
        if (_f_GbbCC->fillnext("CF_Gbb_0l_B", {mBjets_Tmin >= 130*GeV, etmiss > 550*GeV, smallJets[0].pt() > 65*GeV, meff > 2600*GeV})){
          _c["SR_Gbb_0l_B"]->fill();
        }
        else if (mBjets_Tmin >= 130*GeV && etmiss > 550*GeV && smallJets[0].pt() > 65*GeV && meff < 2400*GeV){
          _c["VR_Gbb_0l_B"]->fill();
        }
        if (_f_GbbCC->fillnext("CF_Gbb_0l_M", {mBjets_Tmin >= 130*GeV, etmiss > 550*GeV, meff > 2000*GeV})){
          _c["SR_Gbb_0l_M"]->fill();
        }
        else if (mBjets_Tmin >= 130*GeV && etmiss < 500*GeV && meff > 2000*GeV){
          _c["VR_Gbb_0l_M"]->fill();
        }
        if (_f_GbbCC->fillnext("CF_Gbb_0l_C", {mBjets_Tmin >= 130*GeV, etmiss > 550*GeV, meff > 1600*GeV})){
          _c["SR_Gbb_0l_C"]->fill();
        }
        else if (mBjets_Tmin >= 130*GeV && etmiss < 460*GeV && meff > 1600*GeV){
          _c["VR_Gbb_0l_C"]->fill();
        }
        // CC Gtb 0l SR & VRss
        if (_f_GtbCC->fillnext("CF_Gtb_0l_B", {mBjets_Tmin > 130*GeV, meff>=2500*GeV, etmiss >= 550*GeV, totJetMass >= 200*GeV})){
          _c["SR_Gtb_0l_B"]->fill();
        }
        else if (mBjets_Tmin > 130*GeV && meff<2500*GeV && etmiss >= 550*GeV && totJetMass >= 200*GeV){
          _c["VR_Gtb_0l_B"]->fill();
        }
        if (_f_GtbCC->fillnext("CF_Gtb_0l_M", {mBjets_Tmin > 130*GeV, njets >= 6, nbs >= 4, meff >= 2000*GeV, etmiss >= 550*GeV, totJetMass >= 200*GeV})){
          _c["SR_Gtb_0l_M"]->fill();
        }
        else if (mBjets_Tmin > 130*GeV && njets >= 6 && nbs >= 4 && meff >= 2000*GeV && etmiss < 550*GeV && totJetMass >= 200*GeV){
          _c["VR_Gtb_0l_M"]->fill();
        }
        if (_f_GtbCC->fillnext("CF_Gtb_0l_C", {mBjets_Tmin > 130*GeV, njets >= 7, nbs >= 4, meff >= 1300*GeV, etmiss >= 500*GeV, totJetMass >= 50*GeV})){
          _c["SR_Gtb_0l_C"]->fill();
        }
        else if (mBjets_Tmin > 130*GeV && njets >= 7 && nbs >= 4 && meff >= 1300*GeV && etmiss < 500*GeV && totJetMass >= 50*GeV){
          _c["VR_Gtb_0l_C"]->fill();
        }
      }
      /////////////////////////////////////////////////////////////////////////
      // CC 1l channel
      /////////////////////////////////////////////////////////////////////////
      else if (sigleptons.size() >= 1){
        // 1l channel only variable:
        const double m_transverse = mT(ETMiss.vectorEtMiss(), sigleptons[0].p3());

        // CC Gtt 1l channel
        if (_f_Gtt1lCC->fillnext("CF_Gtt_1l_B", {true, njets>=4 && nbs >= 3, etmiss >= 600*GeV, meff >= 2300*GeV, m_transverse >= 150*GeV, mBjets_Tmin > 120*GeV, totJetMass >= 200*GeV})){
          _c["SR_Gtt_1l_B"]->fill();
        }
        else if (njets>=4 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1500*GeV && m_transverse >= 150*GeV && totJetMass < 200*GeV){
          _c["VR1_Gtt_1l_B"]->fill();
        }
        else if (njets>=5 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1200*GeV && m_transverse < 150*GeV && totJetMass >= 200*GeV){
          _c["VR2_Gtt_1l_B"]->fill();
        }
        else if (njets==4 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1500*GeV && m_transverse < 150*GeV){
          _c["CR_Gtt_1l_B"]->fill();
        }
        
        if (_f_Gtt1lCC->fillnext("CF_Gtt_1l_M1", {true, njets>=5 && nbs >= 3, etmiss >= 600*GeV, meff >= 2000*GeV, m_transverse >= 200*GeV, mBjets_Tmin > 120*GeV, totJetMass >= 200*GeV})){
          _c["SR_Gtt_1l_M1"]->fill();
        }
        else if (njets>=5 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1200*GeV && m_transverse >= 200*GeV && totJetMass < 200*GeV){
          _c["VR1_Gtt_1l_M1"]->fill();
        }
        else if (njets>=6 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1200*GeV && m_transverse < 200*GeV && mBjets_Tmin >= 120*GeV && totJetMass >= 100*GeV){
          _c["VR2_Gtt_1l_M1"]->fill();
        }
        else if (njets==5 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1200*GeV && m_transverse < 200*GeV){
          _c["CR_Gtt_1l_M1"]->fill();
        }

        if (_f_Gtt1lCC->fillnext("CF_Gtt_1l_M2", {true, njets>=8 && nbs >= 3, etmiss >= 500*GeV, meff >= 1100*GeV, m_transverse >= 200*GeV, mBjets_Tmin > 120*GeV, totJetMass >= 100*GeV})){
          _c["SR_Gtt_1l_M2"]->fill();
        }
        else if (njets>=8 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse >= 200*GeV && totJetMass < 100*GeV){
          _c["VR1_Gtt_1l_M2"]->fill();
        }
        else if (njets>=9 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse < 200*GeV && mBjets_Tmin >= 120*GeV && totJetMass < 100*GeV){
          _c["VR2_Gtt_1l_M2"]->fill();
        }
        else if (njets==8 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse < 200*GeV){
          _c["CR_Gtt_1l_M2"]->fill();
        }

        if (_f_Gtt1lCC->fillnext("CF_Gtt_1l_C", {true, njets>=9 && nbs >= 3, etmiss >= 300*GeV, meff >= 800*GeV, m_transverse >= 150*GeV, mBjets_Tmin > 120*GeV})){
          _c["SR_Gtt_1l_C"]->fill();
        }
        else if (njets>=9 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse >= 150*GeV && mBjets_Tmin < 120*GeV){
          _c["VR1_Gtt_1l_C"]->fill();
        }
        else if (njets>=10 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse < 150*GeV && mBjets_Tmin >= 120*GeV){
          _c["VR2_Gtt_1l_C"]->fill();
        }
        else if (njets==9 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse < 150*GeV){
          _c["CR_Gtt_1l_C"]->fill();
        }

        // ///////////////////////////////////////////////////////////////////////
        // There are also some 0l channel control regions with one sig lep:
        // Gtt_0l CRs
        if (njets >= 4 && nbs >= 3 && etmiss >= 200*GeV && meff >= 2000*GeV && m_transverse < 150*GeV && totJetMass >= 150*GeV){
          _c["CR_Gtt_0l_B"]->fill();
        }
        if (njets >= 8 && nbs >= 3 && etmiss >= 200*GeV && meff >= 1100*GeV && m_transverse < 150*GeV && totJetMass >= 150*GeV){
          _c["CR_Gtt_0l_M1"]->fill();
        }
        if (njets >= 9 && nbs >= 3 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse < 150*GeV && totJetMass >= 100*GeV){
          _c["CR_Gtt_0l_M2"]->fill();
        }
        if (njets >= 9 && nbs >= 4 && etmiss >= 200*GeV && meff >= 800*GeV && m_transverse < 150*GeV && totJetMass >= 100*GeV){
          _c["CR_Gtt_0l_C"]->fill();
        }

        // Gbb 0l CRs
        if (smallJets[0].pt() > 65*GeV && meff > 2600*GeV && etmiss > 450*GeV && m_transverse < 150*GeV){
          _c["CR_Gbb_0l_B"]->fill();
        }
        if (meff > 2000*GeV && etmiss > 550*GeV && m_transverse < 150*GeV){
          _c["CR_Gbb_0l_M"]->fill();
        }
        if (meff > 1600*GeV && etmiss > 550*GeV && m_transverse < 150*GeV){
          _c["CR_Gbb_0l_C"]->fill();
        }
        // Gtb 0l CRs
        if (njets >= 4 && nbs >= 3 && meff > 2200*GeV && etmiss > 400*GeV && m_transverse < 150*GeV && totJetMass > 200*GeV){
          _c["CR_Gtb_0l_B"]->fill();
        }
        if (njets >= 6 && nbs >= 4 && meff > 1700*GeV && etmiss > 300*GeV && m_transverse < 150*GeV && totJetMass > 200*GeV){
          _c["CR_Gtb_0l_M"]->fill();
        }
        if (njets >= 7 && nbs >= 4 && meff > 1300*GeV && etmiss > 350*GeV && m_transverse < 150*GeV && totJetMass > 50*GeV){
          _c["CR_Gtb_0l_C"]->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = 1000*crossSection()*(139)/sumW();
      MSG_DEBUG("scale factor " << sf);
      //Rescale counters
      for (const pair<string,CounterPtr> cntr : _c){
        if (cntr.second->sumW() != 0.0){
          cntr.second->scaleW(sf);
        }
        // DEBUG only
        MSG_DEBUG(cntr.first << ": " << cntr.second->sumW() << " +/- " << cntr.second->err() << " (" << cntr.second->numEntries() << " fills)");
      }

      // Rescale cutflows.
      _f_GttNN->scale(sf); _f_GbbNN->scale(sf);
      _f_GtbCC->scale(sf); _f_GbbCC->scale(sf);
      _f_Gtt0lCC->scale(sf); _f_Gtt1lCC->scale(sf);

      // Print cutflows
      MSG_DEBUG("Scaled Flows:\n" << _f_GttNN << _f_GbbNN << _f_GtbCC << _f_GbbCC << _f_Gtt1lCC << _f_Gtt0lCC);
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    /// @}

    private:

    /// @name Useful/conveniece functions for this analysis (static)
    /// @{

    // Normalise the kinematic input (i.e. 1st 84 variables) to the NN
    // TODO: logically, this is static, but I can't quite work out how to square this with the member variables.
    void normalise_kinematic_nn_input(vector<float> &nn_input) const { 
      for (size_t i = 0; i < 84; ++i){ // Magic 84 comes from 84 kinematic variables.
        nn_input[i] = ((nn_input[i] - _nn_norm_means[i])/_nn_norm_devs[i]);
      }
      return;
    }

    // Normalise the paramter nn inputs (last 3 variables)
    void normalise_parameter_nn_input(vector<float> & nn_input) const {
      for (size_t i = 84; i < 87; ++i){ // Magic numbers - last 3 elements of 
        nn_input[i] = ((nn_input[i] - _nn_norm_means[i])/_nn_norm_devs[i]);
      }
      return;
    }

    // normalise just lepton and met inputs (useful for CRZ regions)
    void normalise_lepmet_kinematic_nn_input(vector<float> & nn_input) const {
      for (size_t i = 66; i < 84; ++i){
        nn_input[i] = ((nn_input[i] - _nn_norm_means[i])/_nn_norm_devs[i]);
      }
      return;
    }

    // Zero-out lepton variables (for CRZ regions)
    static void zero_lepton_kinematic_nn_input(vector<float> & nn_input) {
      for (size_t i = 66; i < 82; ++i){
        nn_input[i] = 0.0;
      }
      return;
    }

    // Based on ATL-PHYS-PUB-2016-012 fig 12
    static double JET_BTAG_MV2C10_77_2016(const Jet & j){
      if (j.bTagged()){
        // N.B. !!! There is no overflow bin.
        // There is a clear downward curve in efficiency as pT->inf, but the value to use for pT > 500
        // is super-unclear. The used value (0.69) is certainly a reasonable extrapolation, but it should
        // be noted that tuning this value can make a really significant difference to the bin counts O(10%)
        // HOWEVER: atlas really needs to start increasing the range of this sort of plot
        // The behaviour of jets pT>500 is much more interesting and significant than those in the range 20<pt<25.
        const static vector<double>binedges_pt = {0.00, 30.0, 40.00, 50.00, 60.0, 75.00, 90.0, 105., 150., 200., 500 };
        const static vector<double> eff_pt = {0.63, 0.705, 0.74, 0.76, 0.775, 0.785, 0.795, 0.80, 0.79, 0.75, 0.68};

        return eff_pt[binIndex(j.pt(), binedges_pt, true)];
      }
      else if (j.cTagged()){
        //n.b. there is also a pt/eta table for c mistags, but if that's significant I'll eat my hat.
        return 1./5.;
      }
      else if (j.tauTagged()){
        return 1./16.;
      }
      else return 1./110.;
    }

    /// @name Other private variables belonging to the analysis
    /// @{
    
    // Lists of signal/control/validation regions for Cut'n'Count and Neural Net based analyses.
    const vector<string> _cc_srs = {"SR_Gtt_0l_B", "SR_Gtt_0l_M1", "SR_Gtt_0l_M2", "SR_Gtt_0l_C", "SR_Gbb_0l_B", "SR_Gbb_0l_M", "SR_Gbb_0l_C", "SR_Gtb_0l_B", "SR_Gtb_0l_M", "SR_Gtb_0l_C", "SR_Gtt_1l_B", "SR_Gtt_1l_M1", "SR_Gtt_1l_M2", "SR_Gtt_1l_C"};
    const vector<string> _cc_crs = {"CR_Gtt_0l_B", "CR_Gtt_0l_M1", "CR_Gtt_0l_M2", "CR_Gtt_0l_C", "CR_Gbb_0l_B", "CR_Gbb_0l_M", "CR_Gbb_0l_C", "CR_Gtb_0l_B", "CR_Gtb_0l_M", "CR_Gtb_0l_C", "CR_Gtt_1l_B", "CR_Gtt_1l_M1", "CR_Gtt_1l_M2", "CR_Gtt_1l_C"};
    const vector<string> _cc_vrs = {"VR_Gtt_0l_B", "VR_Gtt_0l_M1", "VR_Gtt_0l_M2", "VR_Gtt_0l_C", "VR_Gbb_0l_B", "VR_Gbb_0l_M", "VR_Gbb_0l_C", "VR_Gtb_0l_B", "VR_Gtb_0l_M", "VR_Gtb_0l_C", "VR1_Gtt_1l_B", "VR1_Gtt_1l_M1", "VR1_Gtt_1l_M2", "VR1_Gtt_1l_C", "VR2_Gtt_1l_B", "VR2_Gtt_1l_M1", "VR2_Gtt_1l_M2", "VR2_Gtt_1l_C"};
    const vector<string> _nn_srs = {"SR_Gtt_2100_1", "SR_Gtt_1800_1", "SR_Gtt_2300_1200", "SR_Gtt_1900_1400", "SR_Gbb_2800_1400", "SR_Gbb_2300_1000", "SR_Gbb_2100_1600", "SR_Gbb_2000_1800"};
    const vector<string> _nn_crs = {"CR_Gtt_2100_1", "CR_Gtt_1800_1", "CR_Gtt_2300_1200", "CR_Gtt_1900_1400", "CR_Gbb_2800_1400", "CR_Gbb_2300_1000", "CR_Gbb_2100_1600", "CR_Gbb_2000_1800", "CRZ_Gbb_2800_1400", "CRZ_Gbb_2300_1000", "CRZ_Gbb_2100_1600", "CRZ_Gbb_2000_1800"};
    const vector<string> _nn_vrs = {"VR_Gtt_2100_1", "VR_Gtt_1800_1", "VR_Gtt_2300_1200", "VR_Gtt_1900_1400", "VR1_Gbb_2800_1400", "VR1_Gbb_2300_1000", "VR1_Gbb_2100_1600", "VR1_Gbb_2000_1800", "VR2_Gbb_2800_1400", "VR2_Gbb_2300_1000", "VR2_Gbb_2100_1600", "VR2_Gbb_2000_1800"};

    // Neural Net cutflows
    CutflowsPtr _f_GbbNN, _f_GttNN;
    // Cut'n'count cut flows
    CutflowsPtr _f_GbbCC, _f_GtbCC, _f_Gtt0lCC, _f_Gtt1lCC;

    // The neural net
    unique_ptr<RivetONNXrt> _nn;

    //Normalisation constants.
    vector<float> _nn_norm_means;
    vector<float> _nn_norm_devs;
    
    const vector<pair<size_t, size_t>> _gtt_nn_params{{1900,1400}, {1800,1}, {2100,1}, {2300,1200}};
    const vector<pair<size_t, size_t>> _gbb_nn_params{{2000,1800}, {2100,1600}, {2300,1000}, {2800,1400}};

    /// @}

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2182381);
}