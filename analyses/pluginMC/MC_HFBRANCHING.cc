// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


    /// @brief MC analysis to compute semi-leptonic branching ratios of heavy-flavour hadrons
    class MC_HFBRANCHING : public Analysis {
    public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(MC_HFBRANCHING);

      /// @name Analysis methods
      /// @{

      const string hadron_id(const int pid) const {
        switch (pid) {
          case PID::B0:           return "B0";
          case PID::BPLUS:        return "BPLUS";
          case PID::B0S:          return "B0S";
          case PID::LAMBDAB:      return "LAMBDAB";
          case PID::D0:           return "D0";
          case PID::DPLUS:        return "DPLUS";
          case PID::DSPLUS:       return "DSPLUS";
          case PID::LAMBDACPLUS:  return "LAMBDACPLUS";
          default:                return "";
        }
      }

      //Semi-leptonic decays with one hadron
      const vector<int> decay_modes_3body(const int pid) const {
        switch (pid) {
          case PID::B0:           return {413,411,10413,10411,20413,415};
          case PID::BPLUS:        return {423,421,10423,10421,20423,425};
          case PID::B0S:          return {433,431,10433,10431,20433,435};
          case PID::LAMBDAB:      return {4122,14122,4124};
          case PID::D0:           return {323,321,10323,325,211,213};
          case PID::DPLUS:        return {313,311,10313,315,111,113,221,331,231};
          case PID::DSPLUS:       return {333,221,331,311,313};
          case PID::LAMBDACPLUS:  return {3122,3212,3214,2112,2114};
          default:                return {};
        }
      }
      //Semi-leptonic decays with two hadrons
      const vector<int> decay_modes_4body(const int pid) const {
        switch (pid) {
          case PID::D0:           return {321,111,211,311};
          case PID::DPLUS:        return {311,111,321,211};
          case PID::DSPLUS:       return {};
          case PID::LAMBDACPLUS:  return {211,211,111,2112};
          default:                return {};
        }
      }

      void fill_Histos(const string &hadron_type, const Particle &p) {

        //********Find decay products of hadron*******//

        vector<int> decay_par; //Vector with direct descendants in hadron decay
        vector<double> child_pt_LAB, child_pt_COM; //Vector with pT of children from hadron decay
        decay_par.clear(), decay_par.resize(0); //vector with hadron's direct descendants
        child_pt_LAB.clear(), child_pt_LAB.resize(0);   //vector with pT of hadron's direct descendants
        child_pt_COM.clear(), child_pt_COM.resize(0);   //vector with pT of hadron's direct descendants in hadron's COM

        for( const Particle& child : p.children()){
          decay_par.push_back(child.abspid()); //Get a list of the direct descendants from the current particle
          child_pt_LAB.push_back(child.pT()/GeV);
          // Reset the boost
          _boost = combine(_boost, _boost.inverse());
          _boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
          Particle temp = p;
          Particle temp_child = child;
          // temp.setMomentum(_boost.transform(temp.momentum())); // to test that the boost works and gives hadron's pt=0
          temp_child.setMomentum(_boost.transform(temp_child.momentum())); //transform child's pt in hadron's COM

          child_pt_COM.push_back(temp_child.pT()/GeV);
        }
        //remove photons from vector to consider QED radiation effects
        decay_par.erase(std::remove(decay_par.begin(), decay_par.end(), 22), decay_par.end());


        //********Compute branching fractions and fill histograms*******//

        vector<int> Modes_3body = decay_modes_3body(p.abspid()); //vector with decay modes of 3-body decays
        vector<int> Modes_4body = decay_modes_4body(p.abspid()); //some 4-body decays considered for c-hadrons
        //Bins to be filled, different for b- and c-hadrons
        int bin_position = -1;
        int last_bin_position = -1;
        bool found_decay_mode = false;

        if (decay_par.size() == 3){  //Semileptonic decays with exactly three decay products
          int lepton_position = -1; //Lepton position in decay_par vector
          for (unsigned int i=0; i < Modes_3body.size(); i++){
            if(contains(decay_par, Modes_3body[i])){

              _h["pt_"+hadron_id(p.abspid())]->fill(p.pT()/GeV); //hadron's pT

              if(hadron_type=="b"){
                bin_position = 3*i; //For b-hadrons there are e, mu and tau decays
                last_bin_position = 3*Modes_3body.size()+1;
              }
              else if(hadron_type=="c"){
                bin_position = 2*i; //For c-hadrons there are e and mu decays only
                last_bin_position = 2*Modes_3body.size()+Modes_4body.size()+1;
              }
              else cout<<"I compute decays of heavy-flavour hadrons, pass b or c"<<endl;

              //electron decays
              if((contains(decay_par, 11) && contains(decay_par, 12))){
                _h[hadron_id(p.abspid())+"_frac_clnu"]->fill(bin_position+1); //branching fraction
                                                                              //Place the e, mu and tau decays for each mode successively
                found_decay_mode = true;

                lepton_position = std::find(decay_par.begin(), decay_par.end(), 11)-decay_par.begin();
                _h[hadron_id(p.abspid())+"_e_pT"] -> fill(child_pt_COM[lepton_position]); //lepton's pT
                _h[hadron_id(p.abspid())+"_lepton_pT_COM"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_LAB"] -> fill(child_pt_LAB[lepton_position]);
              }
              //muon decays
              else if((contains(decay_par, 13) && contains(decay_par, 14))){
                _h[hadron_id(p.abspid())+"_frac_clnu"]->fill(bin_position+2);
                found_decay_mode = true;
                lepton_position = std::find(decay_par.begin(), decay_par.end(), 13)-decay_par.begin();
                _h[hadron_id(p.abspid())+"_mu_pT"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_COM"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_LAB"] -> fill(child_pt_LAB[lepton_position]);
              }
              //tau decays
              else if((contains(decay_par, 15) && contains(decay_par, 16))){
                _h[hadron_id(p.abspid())+"_frac_clnu"]->fill(bin_position+3);
                found_decay_mode = true;
                lepton_position = std::find(decay_par.begin(), decay_par.end(), 15)-decay_par.begin();
                _h[hadron_id(p.abspid())+"_tau_pT"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_COM"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_LAB"] -> fill(child_pt_LAB[lepton_position]);
              }
            }
          }
        }
        else if (decay_par.size() == 4 && hadron_type=="c"){  //Semileptonic decays with exactly four decay products
          int lepton_position = -1;
          for (unsigned int i=0; i < Modes_4body.size(); i++){
            if(contains(decay_par, Modes_4body[i]) && contains(decay_par, Modes_4body[i+1])){

              _h["pt_"+hadron_id(p.abspid())]->fill(p.pT()/GeV); //hadron's pT

              bin_position = 2*Modes_3body.size()+i; //place the 4-body decay modes after the 3-body ones
              last_bin_position = 2*Modes_3body.size()+Modes_4body.size()+1;

              //electron decays
              if((contains(decay_par, 11) && contains(decay_par, 12))){
                _h[hadron_id(p.abspid())+"_frac_clnu"]->fill(bin_position+1);
                found_decay_mode = true;
                lepton_position = std::find(decay_par.begin(), decay_par.end(), 11)-decay_par.begin();
                _h[hadron_id(p.abspid())+"_e_pT"] -> fill(child_pt_COM[lepton_position]); //lepton's pT
                _h[hadron_id(p.abspid())+"_lepton_pT_COM"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_LAB"] -> fill(child_pt_LAB[lepton_position]);
              }
              //muon decays
              else if((contains(decay_par, 13) && contains(decay_par, 14))){
                _h[hadron_id(p.abspid())+"_frac_clnu"]->fill(bin_position+2);
                found_decay_mode = true;
                lepton_position = std::find(decay_par.begin(), decay_par.end(), 13)-decay_par.begin();
                _h[hadron_id(p.abspid())+"_mu_pT"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_COM"] -> fill(child_pt_COM[lepton_position]);
                _h[hadron_id(p.abspid())+"_lepton_pT_LAB"] -> fill(child_pt_LAB[lepton_position]);
              }
            }
          }
        }
        //Fill last bin of branching fractions with decays that don't fall into any category
        if(!found_decay_mode) _h[hadron_id(p.abspid())+"_frac_clnu"]->fill(last_bin_position);
      }

      /// Book histograms and initialise projections before the run
      void init() {

        declare(HeavyHadrons(Cuts::pT > 5.*GeV && Cuts::abseta < 2.5),"HA");

        // histograms
        //Branching ratios
        //Include semi-leptonic decay modes only
        book(_h["B0_frac_clnu"], "BR_B0_clnu", 19, 0.5,19.5);
        book(_h["B0S_frac_clnu"], "BR_B0S_clnu", 19, 0.5,19.5);
        book(_h["BPLUS_frac_clnu"], "BR_BPLUS_clnu", 19, 0.5,19.5);
        book(_h["LAMBDAB_frac_clnu"], "BR_LAMBDAB_clnu", 10, 0.5,10.5);
        book(_h["D0_frac_clnu"], "BR_D0_clnu", 17, 0.5,17.5);
        book(_h["DPLUS_frac_clnu"], "BR_DPLUS_clnu", 23, 0.5,23.5);
        book(_h["DSPLUS_frac_clnu"], "BR_DSPLUS_clnu", 11, 0.5,11.5);
        book(_h["LAMBDACPLUS_frac_clnu"], "BR_LAMBDACPLUS_clnu", 15, 0.5,15.5);

        //Hadron momentum
        book(_h["pt_B0"         ], "B0_pT",          40, 0., 200.);
        book(_h["pt_BPLUS"      ], "BPLUS_pT",       40, 0., 200.);
        book(_h["pt_B0S"        ], "B0S_pT",         40, 0., 200.);
        book(_h["pt_D0"         ], "D0_pT",          40, 0., 200.);
        book(_h["pt_DPLUS"      ], "DPLUS_pT",       40, 0., 200.);
        book(_h["pt_DSPLUS"     ], "DSPLUS_pT",      40, 0., 200.);
        book(_h["pt_LAMBDAB"    ], "LAMBDAB_pT",     40, 0., 200.);
        book(_h["pt_LAMBDACPLUS"], "LAMBDACPLUS_pT", 40, 0., 200.);

        //Chared lepton momentum in LAB
        book(_h["B0_lepton_pT_LAB"         ], "B0_lepton_pT_LAB",          20, 0., 100.);
        book(_h["BPLUS_lepton_pT_LAB"      ], "BPLUS_lepton_pT_LAB",       20, 0., 100.);
        book(_h["B0S_lepton_pT_LAB"        ], "B0S_lepton_pT_LAB",         20, 0., 100.);
        book(_h["D0_lepton_pT_LAB"         ], "D0_lepton_pT_LAB",          20, 0., 100.);
        book(_h["DPLUS_lepton_pT_LAB"      ], "DPLUS_lepton_pT_LAB",       20, 0., 100.);
        book(_h["DSPLUS_lepton_pT_LAB"     ], "DSPLUS_lepton_pT_LAB",      20, 0., 100.);
        book(_h["LAMBDAB_lepton_pT_LAB"    ], "LAMBDAB_lepton_pT_LAB",     20, 0., 100.);
        book(_h["LAMBDACPLUS_lepton_pT_LAB"], "LAMBDACPLUS_lepton_pT_LAB", 20, 0., 100.);


        //Chared lepton momentum in COM
        book(_h["B0_e_pT"], "B0_e_pT", 25, 0., 2.5);
        book(_h["B0_mu_pT"], "B0_mu_pT", 25, 0., 2.5);
        book(_h["B0_tau_pT"], "B0_tau_pT", 25, 0., 2.5);
        book(_h["B0_lepton_pT_COM"], "B0_lepton_pT_COM", 25, 0., 2.5);

        book(_h["B0S_e_pT"], "B0S_e_pT", 25, 0., 2.5);
        book(_h["B0S_mu_pT"], "B0S_mu_pT", 25, 0., 2.5);
        book(_h["B0S_tau_pT"], "B0S_tau_pT", 25, 0., 2.5);
        book(_h["B0S_lepton_pT_COM"], "B0S_lepton_pT_COM", 25, 0., 2.5);

        book(_h["BPLUS_e_pT"], "BPLUS_e_pT", 25, 0., 2.5);
        book(_h["BPLUS_mu_pT"], "BPLUS_mu_pT", 25, 0., 2.5);
        book(_h["BPLUS_tau_pT"], "BPLUS_tau_pT", 25, 0., 2.5);
        book(_h["BPLUS_lepton_pT_COM"], "BPLUS_lepton_pT_COM", 25, 0., 2.5);

        book(_h["LAMBDAB_e_pT"], "LAMBDAB_e_pT", 25, 0., 2.5);
        book(_h["LAMBDAB_mu_pT"], "LAMBDAB_mu_pT", 25, 0., 2.5);
        book(_h["LAMBDAB_tau_pT"], "LAMBDAB_tau_pT", 25, 0., 2.5);
        book(_h["LAMBDAB_lepton_pT_COM"], "LAMBDAB_lepton_pT_COM", 25, 0., 2.5);

        book(_h["D0_e_pT"], "D0_e_pT", 15, 0., 1.5);
        book(_h["D0_mu_pT"], "D0_mu_pT", 15, 0., 1.5);
        book(_h["D0_lepton_pT_COM"], "D0_lepton_pT_COM", 15, 0., 1.5);

        book(_h["DSPLUS_e_pT"], "DSPLUS_e_pT", 15, 0., 1.5);
        book(_h["DSPLUS_mu_pT"], "DSPLUS_mu_pT", 15, 0., 1.5);
        book(_h["DSPLUS_lepton_pT_COM"], "DSPLUS_lepton_pT_COM", 15, 0., 1.5);

        book(_h["DPLUS_e_pT"], "DPLUS_e_pT", 15, 0., 1.5);
        book(_h["DPLUS_mu_pT"], "DPLUS_mu_pT", 15, 0., 1.5);
        book(_h["DPLUS_lepton_pT_COM"], "DPLUS_lepton_pT_COM", 15, 0., 1.5);

        book(_h["LAMBDACPLUS_e_pT"], "LAMBDACPLUS_e_pT", 25, 0., 2.5);
        book(_h["LAMBDACPLUS_mu_pT"], "LAMBDACPLUS_mu_pT", 25, 0., 2.5);
        book(_h["LAMBDACPLUS_lepton_pT_COM"], "LAMBDACPLUS_lepton_pT_COM", 25, 0., 2.5);
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        const HeavyHadrons &ha = apply<HeavyHadrons>(event,"HA");

        if (ha.bHadrons().empty() && ha.cHadrons().empty()) vetoEvent;

        //b hadrons branching ratios
        for (const Particle &hadron : ha.bHadrons()) {
          //Compute branching ratios for listed b-hadrons
          if (hadron_id(hadron.abspid())!="") fill_Histos("b", hadron);
        }
        //c hadrons branching ratios
        for (const Particle &hadron : ha.cHadrons()) {
          if( !hadron.fromBottom()){ //take into account only c-hadrons that don't come from a b-hadron decay
            //Compute branching ratios for listed c-hadrons
            if (hadron_id(hadron.abspid())!="") fill_Histos("c", hadron);
          }
        }
      } // Close Event

      /// Normalise histograms etc., after the run
      void finalize() {
        normalize(_h);
      }

      /// @}


    private:
      /// @name Histograms
      /// @{
      map<string, Histo1DPtr> _h;
      /// @}
      LorentzTransform _boost;

  };


  RIVET_DECLARE_PLUGIN(MC_HFBRANCHING);

}
