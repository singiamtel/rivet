// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  /// @brief Multijet events at 13 TeV
  class CMS_2022_I2170533 : public Analysis {

     public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2022_I2170533);

      /// @name Analysis methods
      ///@{

      /// Book histograms and initialise projections before the run
      void init() {
          FinalState fs;
          FastJets akt(fs, FastJets::ANTIKT, 0.4);
          declare(akt, "antikT");

          vector<double> edges = {0., 150., 170., 180.};
          book(_h_Mult_ptmax200, edges);
          book(_h_Mult_ptmax400, edges);
          book(_h_Mult_ptmax800, edges);
          for (size_t i = 0; i < _h_Mult_ptmax200->numBins(); ++i) {
            book(_h_Mult_ptmax200->bin(i+1), i+1, 1, 1);
            book(_h_Mult_ptmax400->bin(i+1), i+4, 1, 1);
            book(_h_Mult_ptmax800->bin(i+1), i+7, 1, 1);
          }

          // Pt of the first 4 jets
          book(_h_pT1_n1, 10, 1, 1);
          book(_h_pT2_n2, 11, 1, 1);
          book(_h_pT3_n3, 12, 1, 1);
          book(_h_pT4_n4, 13, 1, 1);
      }

      /// Perform the per-event analysis
      void analyze(const Event& event) {
          // Preselection cuts for the jets
          const Jets& jets = apply<FastJets>(event, "antikT").jetsByPt(Cuts::pT > 20 * GeV && Cuts::absrap < 3.2);  // rapidity and pt preselection |y| < 3.2 pT > 20 GeV

          int njet = 0;  // jet counting

          // cuts after preselection
          if (jets.size() < 2) vetoEvent;                                             // dijet cut
          if ((jets[0].pT() < 200.) or (jets[1].pT() < 100.)) vetoEvent;              // pt cut on 2 leading jets
          if ((fabs(jets[0].rap()) > 2.5) or (fabs(jets[1].rap()) > 2.5)) vetoEvent;  // |y| < 2.5 cut on leading jets

          double dphi = deltaPhi(jets[0].phi(), jets[1].phi()) / degree;

          for (const Jet& j : jets) {
              if (j.pT() > 50. && fabs(j.rap()) < 2.5) njet = njet + 1;  // Cuts on the extrajets
          }

          if (njet > 7) njet = 7;  // Last bin in Multiplicity is inclusive

          if (jets[0].pT() > 200 && jets[0].pT() <= 400) _h_Mult_ptmax200->fill(dphi, njet);
          if (jets[0].pT() > 400 && jets[0].pT() <= 800) _h_Mult_ptmax400->fill(dphi, njet);
          if (jets[0].pT() > 800) _h_Mult_ptmax800->fill(dphi, njet);

          // Filling dijet 1,2 pT (The events are already selected with leadin jetpT > 200 GeV and 2nd jet pT > 100 GeV)
          _h_pT1_n1->fill(jets[0].pT());
          _h_pT2_n2->fill(jets[1].pT());

          // Filling extra jets pT (jet3 and jet4)
          if (njet > 2) _h_pT3_n3->fill(jets[2].pT());
          if (njet > 3) _h_pT4_n4->fill(jets[3].pT());
      }

      /// Normalise histograms etc., after the run
      void finalize() {
          const double sf = crossSection() / picobarn / sumW();
          scale(_h_Mult_ptmax200, sf);
          scale(_h_Mult_ptmax400, sf);
          scale(_h_Mult_ptmax800, sf);
          scale({_h_pT1_n1, _h_pT2_n2, _h_pT3_n3, _h_pT4_n4}, sf);
      }

      ///@}

      /// @name Histograms
      ///@{

      Histo1DGroupPtr _h_Mult_ptmax200, _h_Mult_ptmax400, _h_Mult_ptmax800;
      Histo1DPtr _h_pT1_n1, _h_pT2_n2, _h_pT3_n3, _h_pT4_n4;

      ///@}
  };

  RIVET_DECLARE_PLUGIN(CMS_2022_I2170533);

}  // namespace Rivet
