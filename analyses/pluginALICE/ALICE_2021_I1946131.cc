// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {

  /// @brief Prompt strange-meson production at 5 TeV
  class ALICE_2021_I1946131 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1946131);

    void mkAverage(const string& avgname, const vector<string>& estnames) {
      for (auto& b : _e[avgname]->bins()) {
        double wtotal = 0., vtotal = 0., etotal = 0.;
        for (const string& ename : estnames) {
          const auto& est = _e[ename]->binAt(b.xMid());
          if (!_e[ename]->isVisible(est.index()))  continue;
          const double w = 1.0 / sqr( est.relErrAvg() );
          wtotal += w;
          vtotal += est.val() * w;
          etotal += sqr( est.errAvg() * w );
        }
        b.set(vtotal / wtotal, sqrt(etotal) / wtotal);
      }
    }

    void init() {

      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_CENT_PBPB", "V0M","V0M");

      const UnstableParticles ufsD0(Cuts::absrap < 0.5 && Cuts::pT > 0.*GeV && Cuts::abspid == 421);
      declare(ufsD0, "ufsD0");

      const UnstableParticles ufsDplus(Cuts::absrap < 0.5 && Cuts::pT > 2.*GeV && Cuts::abspid == 411);
      declare(ufsDplus, "ufsDplus");

      const UnstableParticles ufsDstar(Cuts::absrap < 0.5 && Cuts::pT > 2.*GeV && Cuts::abspid == 413);
      declare(ufsDstar, "ufsDstar");

      book(_c["sow_pp5TeV"], "_sow_pp5TeV");

      size_t pt_idx = 1, part_idx = 7, beam_idx = 11, avg_idx = 17;
      for (const string& part : vector<string>{ "D0", "Dplus", "Dstar" }) {
        for (const string& cent : vector<string>{ "00-10", "30-50" }) {
          book(_h[part+"Pt_"+cent], pt_idx++, 1, 1);
          if (part != "D0") {
            const string refname = mkAxisCode(part_idx++, 1, 1);
            const YODA::Estimate1D ref = refData(refname);
            string rname(part+"_D0"+cent);
            book(_h["num_"+rname], "_num_"+rname, ref);
            book(_h["den_"+rname], "_den_"+rname, ref);
            book(_e[rname], refname);
          }
          else {
            book(_c["sow_PbPb5TeV_"+cent], "_sow_PbPb5TeV_"+cent);
            book(_e["average"+cent], avg_idx++, 1, 1);
          }
          size_t offset = cent == "00-10"? 0 : 3;
          const string brefname = mkAxisCode(beam_idx + offset, 1, 1);
          const YODA::Estimate1D bref = refData(brefname);
          string rname(part+"PbPb_pp"+cent);
          book(_h["num_"+rname], "_num_"+rname, bref);
          book(_h["den_"+rname], "_den_"+rname, bref);
          book(_e[rname], brefname);
        }
        ++beam_idx;
      }

      const string RAAname = mkAxisCode(19, 1, 1);
      const Estimate1D& RAAref = refData(RAAname);
      book(_h["num_RAAPbPb_pp"], "_num_RAA", RAAref);
      book(_h["den_RAAPbPb_pp"], "_den_RAA", RAAref);
      book(_e["RAAPbPb_pp"], RAAname);

    }

    void analyze(const Event& event) {

      const ParticlePair& beam = beams();
      string CollSystem = "Empty";
      const double NN = 208;

      if (beam.first.pid() == PID::LEAD && beam.second.pid() == PID::LEAD) {
        CollSystem = "PBPB";
        if (fuzzyEquals(sqrtS()/GeV, 5020*NN, 1E-3)) CollSystem += "5TeV";
      }
      if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON) {
        CollSystem = "PP";
        if (fuzzyEquals(sqrtS()/GeV, 5020, 1E-3)) CollSystem += "5TeV";
      }

      const Particles& particlesD0 = apply<UnstableParticles>(event,"ufsD0").particles();
      const Particles& particlesDplus = apply<UnstableParticles>(event,"ufsDplus").particles();
      const Particles& particlesDstar = apply<UnstableParticles>(event,"ufsDstar").particles();

      if (CollSystem == "PP5TeV") {
        _c["sow_pp5TeV"]->fill();

        for (const Particle& p : particlesD0) {
          if (p.fromBottom()) continue;
          for (const string& cent : vector<string>{ "00-10", "30-50" }) {
            _h["den_D0PbPb_pp"+cent]->fill(p.pT()/GeV);
            _h["den_RAAPbPb_pp"]->fill(cent == "00-10"? 5. : 40.);
          }
        }

        for (const Particle& p : particlesDplus) {
          if (p.fromBottom()) continue;
          for (const string& cent : vector<string>{ "00-10", "30-50" }) {
            _h["den_DplusPbPb_pp"+cent]->fill(p.pT()/GeV);
          }
        }

        for (const Particle& p : particlesDstar) {
          if(p.fromBottom()) continue;
          for (const string& cent : vector<string>{ "00-10", "30-50" }) {
            _h["den_DstarPbPb_pp"+cent]->fill(p.pT()/GeV);
          }
        }
      }

      const CentralityProjection& centProj = apply<CentralityProjection>(event, "V0M");

      const double cent_val = centProj();

      if (cent_val >= 50.) vetoEvent;

      if (CollSystem == "PBPB5TeV") {
        string cent("");
        if (cent_val < 10.)  cent = "00-10";
        else if(cent_val >= 30. && cent_val < 50.)  cent = "30-50";

        if (cent == "")  vetoEvent;
        _c["sow_PbPb5TeV_"+cent]->fill();
        for (const Particle& p : particlesD0) {
          if (p.fromBottom()) continue;
          _h["D0Pt_"+cent]->fill(p.pT()/GeV);
          _h["den_Dplus_D0"+cent]->fill(p.pT()/GeV);
          _h["den_Dstar_D0"+cent]->fill(p.pT()/GeV);
          _h["num_D0PbPb_pp"+cent]->fill(p.pT()/GeV);
          _h["num_RAAPbPb_pp"]->fill(cent == "00-10"? 5. : 40.);
        }

        for (const Particle& p : particlesDplus) {
          if (p.fromBottom()) continue;
          _h["DplusPt_"+cent]->fill(p.pT()/GeV);
          _h["num_Dplus_D0"+cent]->fill(p.pT()/GeV);
          _h["num_DplusPbPb_pp"+cent]->fill(p.pT()/GeV);
        }

        for (const Particle& p : particlesDstar) {
          if (p.fromBottom()) continue;
          _h["DstarPt_"+cent]->fill(p.pT()/GeV);
          _h["num_Dstar_D0"+cent]->fill(p.pT()/GeV);
          _h["num_DstarPbPb_pp"+cent]->fill(p.pT()/GeV);
        }
      }
    }

    void finalize() {

      for (auto& item : _h) {
        if (item.first.find("_D0") == string::npos)  continue;
        if (item.first.substr(0, 4) == "den_" && item.first.find("_pp") != string::npos) {
          if (item.first.find("RAA") != string::npos) {
            item.second->bin(1).scaleW( _n["00-10"]/_c["sow_pp5TeV"]->sumW() );
            item.second->bin(3).scaleW( _n["30-50"]/_c["sow_pp5TeV"]->sumW() );
            continue;
          }
          const string cent(item.first.substr(item.first.length()-5));
          scale(item.second, _n[cent]/_c["sow_pp5TeV"]->sumW());
          continue;
        }
        string sfname = "sow_PbPb5TeV_";
        if (item.first.find("00-10") != string::npos)  sfname += "00-10";
        else                                           sfname += "30-50";
        const double sf = (item.first.find("PbPb") != string::npos)? 1.0 : 0.5;
        if (item.first.find("RAA") != string::npos) {
          item.second->bin(1).scaleW( sf/_c["sow_pp5TeV00-10"]->sumW() );
          item.second->bin(3).scaleW( sf/_c["sow_pp5TeV30-50"]->sumW() );
          continue;
        }
        scale(item.second, sf / _c[sfname]->sumW());
      }

      for (auto& item : _e) {
        if (item.first.find("_") == string::npos)  continue;
        divide(_h["num_"+item.first], _h["den_"+item.first], item.second);
      }

      for (const string& cen : vector<string>{ "00-10", "30-50" }) {
        mkAverage("average"+cen, { "D0PbPb_pp"+cen, "DplusPbPb_pp"+cen, "DstarPbPb_pp"+cen });
      }

    }

    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    map<string, Estimate1DPtr> _e;
    map<string, double> _n{ {"00-10", 1572.}, {"30-50", 264.8} };

  };

  RIVET_DECLARE_PLUGIN(ALICE_2021_I1946131);

}
