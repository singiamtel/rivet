#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class ATLAS_2015_I1393758 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2015_I1393758);


    void init() {

      declare(FastJets(FinalState(), JetAlg::ANTIKT, 0.4), "Jets");

      book(forward_kappa3, 1, 1, 1);
      book(forwardRMS_kappa3, "d02-x01-y01");

      book(central_kappa3, 3, 1, 1);
      book(centralRMS_kappa3, "d04-x01-y01");

      book(forward_kappa5, 5, 1, 1);
      book(forwardRMS_kappa5, "d06-x01-y01");

      book(central_kappa5, 7, 1, 1);
      book(centralRMS_kappa5, "d08-x01-y01");

      book(forward_kappa7, 9, 1, 1);
      book(forwardRMS_kappa7, "d10-x01-y01");

      book(central_kappa7, 11, 1, 1);
      book(centralRMS_kappa7, "d12-x01-y01");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets m_goodJets = apply<JetFinder>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.1);

      if (m_goodJets.size() < 2)        vetoEvent;
      if (m_goodJets[0].pT() < 50*GeV)  vetoEvent;
      if (m_goodJets[1].pT() < 50*GeV)  vetoEvent;
      if (fabs(1.0 - m_goodJets[0].pT()/m_goodJets[1].pT()) > 0.5)  vetoEvent;

      bool check = m_goodJets[0].abseta() < m_goodJets[1].abseta();
      int pos_f = int(check);
      int pos_c = int(!check);

      double kappa3_f   = CalculateJetCharge(m_goodJets[pos_f], 0.3, 0.5, 1.8);
      double kappa5_f   = CalculateJetCharge(m_goodJets[pos_f], 0.5, 0.5, 1.2);
      double kappa7_f   = CalculateJetCharge(m_goodJets[pos_f], 0.7, 0.5, 0.9);
      double pT_f = m_goodJets[pos_f].pT();

      double kappa3_c = CalculateJetCharge(m_goodJets[pos_c], 0.3, 0.5, 1.8);
      double kappa5_c   = CalculateJetCharge(m_goodJets[pos_c], 0.5, 0.5, 1.2);
      double kappa7_c   = CalculateJetCharge(m_goodJets[pos_c], 0.7, 0.5, 0.9);
      double pT_c = m_goodJets[pos_c].pT();

      forward_kappa3->fill(pT_f, kappa3_f);
      forward_kappa5->fill(pT_f, kappa5_f);
      forward_kappa7->fill(pT_f, kappa7_f);

      central_kappa3->fill(pT_c, kappa3_c);
      central_kappa5->fill(pT_c, kappa5_c);
      central_kappa7->fill(pT_c, kappa7_c);
    }

    double CalculateJetCharge(Jet& jet, double kappa=0.5, double pTcut=0.5, double Qmax=1.2) {
      double PTkap = pow(jet.momentum().pT(),kappa);
      double jetcharge = 0.;
      for (const Particle& p : jet.particles()) {
        if (p.pT() < pTcut)  continue;
        if (p.charge3()) jetcharge += (p.charge3()/3.)*pow(p.pT(),kappa)/PTkap;
      }
      //Overflow and underflow
      if (jetcharge > Qmax) jetcharge = Qmax*0.9999;
      if (jetcharge < -Qmax) jetcharge = -Qmax*0.9999;
      return jetcharge;
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      if (numEvents() > 2) {
        for (unsigned int i = 1; i < forward_kappa3->numBins()+1; ++i) {
          double stdv_fkappa3 = forward_kappa3->bin(i).effNumEntries() > 1? forward_kappa3->bin(i).yStdDev() : 0.0;
          //See Eq. 3 for the factor of two: https://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf
          double yerr_fkappa3  = safediv(sqrt(forward_kappa3->bin(i).sumW2()), 2.*forward_kappa3->bin(i).sumW());
          forwardRMS_kappa3->bin(i).set(stdv_fkappa3, yerr_fkappa3);

          double stdv_fkappa5 = forward_kappa5->bin(i).effNumEntries() > 1? forward_kappa5->bin(i).yStdDev() : 0.0;
          double yerr_fkappa5  = safediv(sqrt(forward_kappa5->bin(i).sumW2()), 2.*forward_kappa5->bin(i).sumW());
          forwardRMS_kappa5->bin(i).set(stdv_fkappa5, yerr_fkappa5);

          double stdv_fkappa7 = forward_kappa7->bin(i).effNumEntries() > 1? forward_kappa7->bin(i).yStdDev() : 0.0;
          double yerr_fkappa7  = safediv(sqrt(forward_kappa7->bin(i).sumW2()), 2.*forward_kappa7->bin(i).sumW());
          forwardRMS_kappa7->bin(i).set(stdv_fkappa7, yerr_fkappa7);

          double stdv_ckappa3 = central_kappa3->bin(i).effNumEntries() > 1? central_kappa3->bin(i).yStdDev() : 0.0;
          double yerr_ckappa3  = safediv(sqrt(central_kappa3->bin(i).sumW2()), 2.*central_kappa3->bin(i).sumW());
          centralRMS_kappa3->bin(i).set(stdv_ckappa3, yerr_ckappa3);

          double stdv_ckappa5 = central_kappa5->bin(i).effNumEntries() > 1? central_kappa5->bin(i).yStdDev() : 0.0;
          double yerr_ckappa5  = safediv(sqrt(central_kappa5->bin(i).sumW2()), 2.*central_kappa5->bin(i).sumW());
          centralRMS_kappa5->bin(i).set(stdv_ckappa5, yerr_ckappa5);

          double stdv_ckappa7 = central_kappa7->bin(i).effNumEntries() > 1? central_kappa7->bin(i).yStdDev() : 0.0;
          double yerr_ckappa7  = safediv(sqrt(central_kappa7->bin(i).sumW2()), 2.*central_kappa7->bin(i).sumW());
          centralRMS_kappa7->bin(i).set(stdv_ckappa7, yerr_ckappa7);

        }
      }

    }


  private:

    Profile1DPtr forward_kappa3;
    Profile1DPtr forward_kappa5;
    Profile1DPtr forward_kappa7;

    Profile1DPtr central_kappa3;
    Profile1DPtr central_kappa5;
    Profile1DPtr central_kappa7;

    Estimate1DPtr forwardRMS_kappa3;
    Estimate1DPtr forwardRMS_kappa5;
    Estimate1DPtr forwardRMS_kappa7;

    Estimate1DPtr centralRMS_kappa3;
    Estimate1DPtr centralRMS_kappa5;
    Estimate1DPtr centralRMS_kappa7;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2015_I1393758);

}
