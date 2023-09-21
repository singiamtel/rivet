// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  namespace {

    inline double calcz(const Jet& j, const Particle& p) {
      const double num = j.p3().dot(p.p3());
      const double den = j.p3().mod2();
      return num/den;
    }

    inline double calcptrel(const Jet& j, const Particle& p) {
      const double num = j.p3().cross(p.p3()).mod();
      const double den = j.p3().mod();
      return num/den;
    }

    inline double calcr(const Jet& j, const Particle& p) {
      return deltaR(j.rapidity(), j.phi(), p.rapidity(), p.phi());
    }

    // For annulus area kludge
    /// @todo Improve somehow... need normalisation *without* bin width factors!
    inline double calcrweight(const Jet& j, const Particle& p) {
      size_t nBins_r = 26;
      double bins_r[] = { 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                          0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30,
                          0.35, 0.40, 0.45, 0.50, 0.55, 0.60 };
      double r = calcr(j,p);
      for (size_t bin = 0 ; bin < nBins_r ; bin++) {
        if (r < bins_r[bin+1]) {
          double up = bins_r[bin+1];
          double down = bins_r[bin];
          return ((up-down)/(M_PI*(up*up-down*down)));
        }
      }
      return 1.0;
    }
  }


  class ATLAS_2011_I919017 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I919017);

    /// @}


  public:

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      ChargedFinalState cfs((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  0.3*GeV));
      FastJets trkjets04(cfs, FastJets::ANTIKT, 0.4);
      FastJets trkjets06(cfs, FastJets::ANTIKT, 0.6);
      declare(trkjets04, "Jets04");
      declare(trkjets06, "Jets06");

      // Book histograms
      const vector<string> rapbins{"00_05_", "05_10_", "10_15_", "15_19_"};
      const vector<string> ptbins{"04_06", "06_10", "10_15", "15_24", "24_40"};
      for (const string R : {"04_", "06_"}) {
        size_t iy = (R == "04_")? 1 : 2;
        book(_c[R+"sumw"], "/TMP/"+R+"sumw");
        size_t irap = 0;
        for (const string& rapbin : rapbins) {
          ++irap;
          book(_hinc[R+rapbin], irap, 1, iy);
          size_t ipt = 0;
          for (const string& ptbin : ptbins) {
            ++ipt;
            const string suff = R + rapbin + ptbin;
            if (irap == 1) {
              size_t offset = rapbins.size();
              book(_h["N"+R+"00_19_"+ptbin], offset+ipt, 1, iy);
              offset += ptbins.size()*ptbins.size();
              book(_h["z"+R+"00_19_"+ptbin], offset+ipt, 1, iy);
              offset += ptbins.size()*ptbins.size();
              book(_h["ptrel"+R+"00_19_"+ptbin], offset+ipt, 1, iy);
              offset += ptbins.size()*ptbins.size();
              book(_h["rdA"+R+"00_19_"+ptbin], offset+ipt, 1, iy);
              book(_c["numjets"+R+"00_19_"+ptbin], "/TMP/numjets"+R+"00_19_"+ptbin);
            }
            size_t offset = rapbins.size() + irap*ptbins.size();
            book(_h["N"+suff], offset+ipt, 1, iy);
            offset = rapbins.size() + (irap+ptbins.size())*ptbins.size();
            book(_h["z"+suff], offset+ipt, 1, iy);
            offset = rapbins.size() + (irap+2*ptbins.size())*ptbins.size();
            book(_h["ptrel"+suff], offset+ipt, 1, iy);
            offset = rapbins.size() + (irap+3*ptbins.size())*ptbins.size();
            book(_h["rdA"+suff], offset+ipt, 1, iy);
            book(_c["numjets"+suff], "/TMP/numjets"+suff);
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      doCollection("04", event);
      doCollection("06", event);
    }

    void doCollection(const string& Rval, const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "Jets"+Rval).jets();
      if (jets.empty())  return;
      const string R = Rval+"_";
      _c[R+"sumw"]->fill();
      for (const Jet& j : jets) {
        const double jetpt = j.pT()/GeV;
        const double jetrap = j.absrap();
        string rapbin, ptbin;
        if (jetrap >= 1.9)  continue;

        if (jetrap < 0.5)       rapbin = "00_05_";
        else if (jetrap < 1.0)  rapbin = "05_10_";
        else if (jetrap < 1.5)  rapbin = "10_15_";
        else                    rapbin = "15_19_";
        _hinc[R+rapbin]->fill(jetpt);

        if (inRange(jetpt, 4., 6.))         ptbin = "04_06";
        else if (inRange(jetpt, 6., 10.))   ptbin = "06_10";
        else if (inRange(jetpt, 10., 15.))  ptbin = "10_15";
        else if (inRange(jetpt, 15., 24.))  ptbin = "15_24";
        else if (inRange(jetpt, 24., 40.))  ptbin = "24_40";
        else  continue;
        fillgroup(j, R+"00_19_"+ptbin);
        fillgroup(j, R+rapbin+ptbin);
      } // each jet
    }

    void fillgroup(const Jet& j, const string& suff) {
      _c["numjets"+suff]->fill();
      _h["N"+suff]->fill(j.particles().size());
      for (const Particle& p : j.particles()) {
        _h["z"+suff]->fill(calcz(j, p));
        _h["ptrel"+suff]->fill(calcptrel(j, p));
        _h["rdA"+suff]->fill(calcr(j, p), calcrweight(j, p));
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double xsec = crossSection()/microbarn;

      for (auto& item : _hinc) {
        if (item.first.substr(0,2) == "04") {
          double sf = _c["04_sumw"]->val();
          if (item.first.substr(3,2) == "15")  sf *= 0.8;
          if (sf)  scale(item.second, xsec/sf);
          else     normalize(item.second, 0.0);
        }
        else {
          double sf = _c["06_sumw"]->val();
          if (item.first.substr(3,2) == "15")  sf *= 0.8;
          if (sf)  scale(item.second, xsec/sf);
          else     normalize(item.second, 0.0);
        }
      }

      for (auto& item : _h) {
        const string suff = item.first.substr(item.first.size()-14);
        const double sf = _c["numjets"+suff]->val();
        if (sf)  scale(item.second, 1.0/sf);
        else     normalize(item.second, 0.0);
      }

    }

    /// @}


  private:

    /// @name Histograms
    /// @{

    map<string,CounterPtr> _c;
    map<string,Histo1DPtr> _h, _hinc;

    /// @}

  };


  // This global object acts as a hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2011_I919017);

}
