// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS Azimuthal corellations at 13 TeV
  class CMS_2018_I1643640 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1643640);

    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.4);
      declare(akt, "antikT");

      const vector<double> edges_2J{200., 300., 400., 500., 600., 700., 800., 1000., 1200., 7000.};
      const vector<double> edges{200., 300., 400., 500., 600., 700., 800., 1000., 7000.};
      book(_h_deltaPhi_2J_phi12, edges_2J);
      book(_h_deltaPhi_3J_phi12, edges);
      book(_h_deltaPhi_4J_phi12, edges);
      book(_h_deltaPhi_3J_phimin2J, edges);
      book(_h_deltaPhi_4J_phimin2J, edges);
      for (size_t i=1; i<_h_deltaPhi_2J_phi12->numBins()+1; ++i) {
        book(_h_deltaPhi_2J_phi12->bin(i), i, 1, 1);
        if (i==9)  continue;
        book(_h_deltaPhi_3J_phi12->bin(i), 9+i, 1, 1);
        book(_h_deltaPhi_4J_phi12->bin(i), 17+i, 1, 1);
        book(_h_deltaPhi_3J_phimin2J->bin(i), 25+i, 1, 1);
        book(_h_deltaPhi_4J_phimin2J->bin(i), 33+i, 1, 1);
      }
    }


    void analyze(const Event & event) {
      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt();

      // 2 jet case and Delta_phi12
      if( jets.size() >= 2 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV) ) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) ) {
            double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
            _h_deltaPhi_2J_phi12->fill(jets[0].pT(), dphi);
          }
        }
      }

      // 3 jet case and Delta_phi12
      if ( jets.size() >= 3 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV) ) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5)) {
            double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
            _h_deltaPhi_3J_phi12->fill(jets[0].pT(), dphi);
          }
        }
      }

      // 4 jet case and Delta_phi12
      if ( jets.size() >= 4 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV)   && (jets[3].pT() >= 100.*GeV)) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5) &&  (fabs(jets[3].rap()) <= 2.5)) {
            double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
            _h_deltaPhi_4J_phi12->fill(jets[0].pT(), dphi);
          }
        }
      }

      // 3 jet case and Delta_Phi_min2j
      if ( jets.size() >= 3 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV) ) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5)) {
            double dphi01 = deltaPhi(jets[0].phi(), jets[1].phi());
            if (dphi01 >= PI/2. ){
              double dphi02 = deltaPhi(jets[0].phi(), jets[2].phi());
              double dphi12 = deltaPhi(jets[1].phi(), jets[2].phi());
              // evaluate DPhi2Jmin
              vector<double> Dphis2J{dphi01,dphi02,dphi12};
              double DPhi2Jmin = min(Dphis2J);
              // double Dphis2J[3] = {dphi01,dphi02,dphi12};
              // double DPhi2Jmin = Dphis2J[0];
              // for (int gg=1; gg<3; ++gg) { if (DPhi2Jmin>Dphis2J[gg]) DPhi2Jmin = Dphis2J[gg]; }
              _h_deltaPhi_3J_phimin2J->fill(jets[0].pT(), DPhi2Jmin);
            }
          }
        }
      }

      // 4 jet case and Delta_Phi_min2j
      if ( jets.size() >= 4 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV)   && (jets[3].pT() >= 100.*GeV)) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5) &&  (fabs(jets[3].rap()) <= 2.5)) {
            double dphi01 = deltaPhi(jets[0].phi(), jets[1].phi());
            if (dphi01 >= PI/2.) {
              double dphi02 = deltaPhi(jets[0].phi(), jets[2].phi());
              double dphi03 = deltaPhi(jets[0].phi(), jets[3].phi());
              double dphi12 = deltaPhi(jets[1].phi(), jets[2].phi());
              double dphi13 = deltaPhi(jets[1].phi(), jets[3].phi());
              double dphi23 = deltaPhi(jets[2].phi(), jets[3].phi());
              /// evaluate DPhi2Jmin
              // double Dphis2J[6]={dphi01,dphi02,dphi03,dphi12,dphi13,dphi23};
              // double DPhi2Jmin=Dphis2J[0];
              // for(int gg=1; gg<6; ++gg){ if(DPhi2Jmin>Dphis2J[gg]){DPhi2Jmin=Dphis2J[gg];} }
              vector<double> Dphis2J{dphi01,dphi02,dphi03,dphi12,dphi13,dphi23};
              double DPhi2Jmin = min(Dphis2J);
              _h_deltaPhi_4J_phimin2J->fill(jets[0].pT(), DPhi2Jmin);
            }
          }
        }
      }
    }  // end analyze


    void finalize() {
      normalize(_h_deltaPhi_2J_phi12);
      normalize(_h_deltaPhi_3J_phi12);
      normalize(_h_deltaPhi_4J_phi12);
      normalize(_h_deltaPhi_3J_phimin2J);
      normalize(_h_deltaPhi_4J_phimin2J);
    }


  private:

    Histo1DGroupPtr _h_deltaPhi_2J_phi12;
    Histo1DGroupPtr _h_deltaPhi_3J_phi12;
    Histo1DGroupPtr _h_deltaPhi_4J_phi12;
    Histo1DGroupPtr _h_deltaPhi_3J_phimin2J;
    Histo1DGroupPtr _h_deltaPhi_4J_phimin2J;

  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1643640);

}
