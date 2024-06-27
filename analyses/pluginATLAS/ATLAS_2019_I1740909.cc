#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief jet fragmentation at 13 TeV
  class ATLAS_2019_I1740909: public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1740909);

    /// Book cuts and projections
    void init() {

      const FinalState bare_MU(Cuts::abspid == PID::MUON);

      VetoedFinalState jetinput;
      jetinput.addVetoOnThisFinalState(bare_MU);
      jetinput.addVetoOnThisFinalState(InvisibleFinalState());

      FastJets jetpro(jetinput, JetAlg::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      book(_p["nch_jetpt_F"], 1, 1, 1);
      book(_p["nch_jetpt_C"], 2, 1, 1);
      book(_p["nch_jetpt_B"], 9, 1, 1);

      vector<double> ptbins{100., 200., 300., 400., 500., 600., 700., 800.,
                            900., 1000., 1200., 1400., 1600., 2000., 2500.};
      for (const string& suff : vector<string>{"B", "F", "C"}) {
        book(_h["nch_"+suff], ptbins);
        book(_h["r_"+suff], ptbins);
        book(_h["zeta_"+suff], ptbins);
        book(_h["pTrel_"+suff], ptbins);
      }
      for (size_t idx = 0; idx < ptbins.size()-1; ++idx) {
        dualbook(_h["nch_B"]->bin(idx+1),    13 + idx, 1, 1);
        dualbook(_h["r_B"]->bin(idx+1),      27 + idx, 1, 1);
        dualbook(_h["zeta_B"]->bin(idx+1),   41 + idx, 1, 1);
        dualbook(_h["pTrel_B"]->bin(idx+1),  55 + idx, 1, 1);
        dualbook(_h["nch_F"]->bin(idx+1),    69 + idx, 1, 1);
        dualbook(_h["r_F"]->bin(idx+1),      83 + idx, 1, 1);
        dualbook(_h["zeta_F"]->bin(idx+1),   97 + idx, 1, 1);
        dualbook(_h["pTrel_F"]->bin(idx+1), 111 + idx, 1, 1);
        dualbook(_h["nch_C"]->bin(idx+1),   125 + idx, 1, 1);
        dualbook(_h["r_C"]->bin(idx+1),     139 + idx, 1, 1);
        dualbook(_h["zeta_C"]->bin(idx+1),  153 + idx, 1, 1);
        dualbook(_h["pTrel_C"]->bin(idx+1), 167 + idx, 1, 1);
      }
    }

    void dualbook(Histo1DPtr& hist, unsigned int d, unsigned int x, unsigned int y) {
      const string hname = mkAxisCode(d, x, y);
      book(hist, "_aux_"+hname, refData(hname));
      book(_s[hist.get()->basePath()], hname);
    }

    void analyze(const Event& event) {

      //Init
      double fnch=0;
      double cnch=0;
      double fzval=0;
      double czval=0;
      double frval=0;
      double crval=0;
      double ftval=0;
      double ctval=0;

      // Event selection
      Jets m_goodJets = apply<JetFinder>(event, "Jets").jetsByPt(Cuts::pT > 100*GeV && Cuts::abseta < 2.1);
      if (m_goodJets.size() < 2) vetoEvent;
      if (fabs(1.0 - m_goodJets[0].pT()/m_goodJets[1].pT()) > 0.5)  vetoEvent;
      // Decide forward or central
      bool check = m_goodJets[0].abseta() < m_goodJets[1].abseta();
      int pos_f = int(check);
      int pos_c = int(!check);

      // Calculate obs, separately for central and fwd, also get bin
      double fpt = m_goodJets[pos_f].pT();
      double cpt = m_goodJets[pos_c].pT();

      for (const Particle& p : m_goodJets[pos_f].particles()) {
        if (p.pT() < 0.5*GeV)  continue;
        if (p.charge() != 0) {
          ++fnch;

          fzval = p.pT() / m_goodJets[pos_f].pt();
          ftval = p.pT()*sin(p.phi()-m_goodJets[pos_f].phi());
          frval = deltaR(m_goodJets[pos_f], p);

          _h["r_F"]->fill(fpt, frval);
          _h["r_B"]->fill(fpt, frval);
          _h["zeta_F"]->fill(fpt, fzval);
          _h["zeta_B"]->fill(fpt, fzval);
          _h["pTrel_F"]->fill(fpt, ftval);
          _h["pTrel_B"]->fill(fpt, ftval);
        }
      }


	    for (const Particle& p : m_goodJets[pos_c].particles()) {
        if (p.pT() < 0.5*GeV)  continue;
        if (p.charge() != 0) {
          ++cnch;
          czval = p.pT() / m_goodJets[pos_c].pt();
          ctval = p.pT()*sin(p.phi()-m_goodJets[pos_c].phi());
          crval = deltaR(m_goodJets[pos_c], p);

          _h["r_C"]->fill(cpt, crval);
          _h["r_B"]->fill(cpt, crval);
          _h["zeta_C"]->fill(cpt, czval);
          _h["zeta_B"]->fill(cpt, czval);
          _h["pTrel_C"]->fill(cpt, ctval);
          _h["pTrel_B"]->fill(cpt, ctval);
        }
      }

      if (fnch > 63)  fnch = 63;
      if (cnch > 63)  cnch = 63;

       //Fill nchg histo

      _p["nch_jetpt_F"]->fill(fpt,fnch);
      _p["nch_jetpt_C"]->fill(cpt,cnch);
      _p["nch_jetpt_B"]->fill(fpt,fnch);
      _p["nch_jetpt_B"]->fill(cpt,cnch);

      _h["nch_F"]->fill(fpt, fnch);
      _h["nch_C"]->fill(cpt, cnch);
      _h["nch_B"]->fill(fpt, fnch);
      _h["nch_B"]->fill(cpt, cnch);
    }


    void finalize() {

      for (const string& suff : vector<string>{"B", "F", "C"}) {
        const double num = suff=="B"? 2.0 : 1.0;
        vector<double> sf = _h["nch_"+suff]->sumWGroup();
        for (double& f : sf) {
          f = safediv(num, f, 0.0);
        }
        scale({_h["nch_"+suff], _h["zeta_"+suff], _h["pTrel_"+suff], _h["r_"+suff]}, sf);
        for (auto& r : _h["r_"+suff]->bins()) {
          for (auto& b : r->bins()) {
            b.scaleW(1.0/(2*M_PI*b.xMid()));
          }
        }
      }

      for (auto& hist : _h) {
        for (auto& b : hist.second->bins()) {
          barchart(b, _s[b.get()->basePath()]);
        }
      }

    }

  private:

    map<string, Histo1DGroupPtr> _h;
    map<string, Estimate1DPtr> _s;
    map<string, Profile1DPtr> _p;
  };

  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1740909);

}
