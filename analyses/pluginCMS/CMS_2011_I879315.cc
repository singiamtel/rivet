// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/HistoGroup.hh"
using namespace std;

namespace Rivet {


  /// Measurement of the NSD charged particle multiplicity at 0.9, 2.36, and 7 TeV
  class CMS_2011_I879315 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I879315);


    void init() {
      ChargedFinalState cfs(Cuts::abseta < 2.4);
      declare(cfs, "CFS");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        size_t offset = 0;
        if (en == "2360"s)  offset = 5;
        else if (en == "7000"s)  offset = 10;
        book(_g[en+"dNch_dn"], {0.0, 0.5, 1.0, 1.5, 2.0, 2.4});
        for (auto& b : _g[en+"dNch_dn"]->bins()) {
          book(b, b.index()+offset+1, 1, 1);
        }

        if (en == "2360"s)  offset = 1;
        else if (en == "7000"s)  offset = 2;
        book(_h[en+"dNch_dn_pt500_eta24"], 20+offset, 1, 1);
        book(_p[en+"dmpt_dNch_eta24"],     23+offset, 1, 1);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    void analyze(const Event& event) {

      // Get the charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");

      // Resetting the multiplicity for the event to 0;
      size_t nBins = _g[_sqs+"dNch_dn"]->numBins();
      vector<int> _nch_in_Evt(nBins, 0);
      vector<int> _nch_in_Evt_pt500(nBins, 0);
      double sumpt = 0;

      // Loop over particles in event
      for (const Particle& p : charged.particles()) {
        // Selecting only charged hadrons
        if (! PID::isHadron(p.pid())) continue;

        const double pT = p.pT(); sumpt += pT;
        const double abseta = p.abseta();
        size_t ieta = nBins;
        while (ieta--) {
          if (abseta > _g[_sqs+"dNch_dn"]->bin(ieta+1).xMax()) break;
          ++_nch_in_Evt[ieta];
          if (pT > 0.5*GeV) ++_nch_in_Evt_pt500[ieta];
        }
      }

      // Filling multiplicity-dependent histogramms
      for (auto& b : _g[_sqs+"dNch_dn"]->bins()) {
        b->fill(_nch_in_Evt[b.index()-1]);
      }

      if (_nch_in_Evt[nBins-1] != 0) {
        _p[_sqs+"dmpt_dNch_eta24"]->fill(_nch_in_Evt[nBins-1], sumpt/GeV / _nch_in_Evt[nBins-1]);
      }
      _h[_sqs+"dNch_dn_pt500_eta24"]->fill(_nch_in_Evt_pt500[nBins-1]);
    }


    void finalize() {
      normalize(_g);
      normalize(_h);
    }


  private:

    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, Histo1DGroupPtr> _g;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_I879315, CMS_2011_S8884919);

}
