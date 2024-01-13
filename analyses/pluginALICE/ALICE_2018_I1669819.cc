// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Strange meson production at 5 TeV
  class ALICE_2018_I1669819 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2018_I1669819);

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

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      const ALICE::PrimaryParticles app(Cuts::abseta < 0.8 && Cuts::pT > 1*GeV && Cuts::abscharge > 0);
      declare(app, "app");

      const UnstableParticles ufsD0(Cuts::absrap < 0.5 && Cuts::pT > 1*GeV && Cuts::abspid == PID::D0);
      declare(ufsD0, "ufsD0");

      const UnstableParticles ufsDplus(Cuts::absrap < 0.5 && Cuts::pT > 2*GeV && Cuts::abspid == PID::DPLUS);
      declare(ufsDplus, "ufsDplus");

      const UnstableParticles ufsDstar(Cuts::absrap < 0.5 && Cuts::pT > 1*GeV && Cuts::abspid == PID::DSTARPLUS);
      declare(ufsDstar, "ufsDstar");

      const UnstableParticles ufsDs(Cuts::absrap < 0.5 && Cuts::pT > 2*GeV && Cuts::abspid == PID::DSPLUS);
      declare(ufsDs, "ufsDs");

      book(_c["sow_pp5TeV"], "_sow_pp5TeV");
      book(_c["sow_pp2TeV"], "_sow_pp2TeV");
      book(_c["sow_PbPb2TeV"], "_sow_PbPb2TeV");

      unsigned int pt_idx = 1, part_idx = 13, beam_idx = 25;
      for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
        book(_c["sow_PbPb5TeV_"+cen], "_sow_PbPb5TeV_"+cen);
        for (const string& part : vector<string>{ "D0", "Dplus", "Dstar", "Ds" }) {

          book(_h[part+"Pt_"+cen], pt_idx++, 1, 1);

          if (part != "D0") {
            const string ref1name = mkAxisCode(part_idx++, 1, 1);
            const YODA::Estimate1D ref1 = refData(ref1name);
            string rname(part+"_D0"+cen);
            book(_h["num_"+rname], "_num_"+rname, ref1);
            book(_h["den_"+rname], "_den_"+rname, ref1);
            book(_e[rname], ref1name);
            if (part == "Ds") {
              const string ref2name = mkAxisCode(part_idx++, 1, 1);
              const YODA::Estimate1D ref2 = refData(ref2name);
              rname = part+"_Dplus_"+cen;
              book(_h["num_"+rname], "_num_"+rname, ref2);
              book(_h["den_"+rname], "_den_"+rname, ref2);
              book(_e[rname], ref2name);
            }
          }
          const string brefname = mkAxisCode(beam_idx++, 1, 1);
          const YODA::Estimate1D bref = refData(brefname);
          string rname(part+"PbPb_pp"+cen);
          book(_h["num_"+rname], "_num_"+rname, bref);
          book(_h["den_"+rname], "_den_"+rname, bref);
          book(_e[rname], brefname);

          if (part == "Ds") {
            book(_e["average"+cen], beam_idx++, 1, 1);
          }
        }
      }

      //Table 40 is PbPb 2.76 TeV
      for (const string& part : vector<string>{ "D0", "Dplus", "Dstar" }) {
        vector<double> bins;
        if (part != "D0") bins = {3., 4., 5., 6., 8., 12., 16., 24., 36.};
        else              bins = {1., 2., 3., 4., 5., 6., 8., 12., 16., 24.};
        book(_h["num_"+part+"PbPb_pp2TeV"], "_"+part+"_PbPb", bins);
        book(_h["den_"+part+"PbPb_pp2TeV"], "_"+part+"_pp",   bins);
        book(_e[part+"PbPb_pp2TeV"], "_"+part+"_PbPb_pp", bins);
      }
      book(_e["average2TeV"], 40, 1, 1);

      unsigned int avg_idx = 41;
      for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
        const string refname = mkAxisCode(avg_idx++, 1, 1);
        const Estimate1D& ref = refData(refname);
        book(_h["num_mult_PbPb_pp"+cen], "_num_mult_PbPb_pp_"+cen, ref);
        book(_h["den_mult_PbPb_pp"+cen], "_den_mult_PbPb_pp_"+cen, ref);
        book(_e["mult_PbPb_pp"+cen], "_mult_PbPb_pp_"+cen);
        book(_e["avgDmult"+cen], refname);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ParticlePair& beam = beams();
      string CollSystem = "Empty";
      const double NN = 208;

      if (beam.first.pid() == PID::LEAD && beam.second.pid() == PID::LEAD) {
        CollSystem = "PBPB";
        if(fuzzyEquals(sqrtS()/GeV, 2760*NN, 1E-3)) CollSystem += "2TeV";
        else if(fuzzyEquals(sqrtS()/GeV, 5020*NN, 1E-3)) CollSystem += "5TeV";
      }
      if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON) {
        CollSystem = "PP";
        if(fuzzyEquals(sqrtS()/GeV, 2760, 1E-3)) CollSystem += "2TeV";
        else if(fuzzyEquals(sqrtS()/GeV, 5020, 1E-3)) CollSystem += "5TeV";
      }

      const Particles& particlesD0 = apply<UnstableParticles>(event,"ufsD0").particles();
      const Particles& particlesDplus = apply<UnstableParticles>(event,"ufsDplus").particles();
      const Particles& particlesDstar = apply<UnstableParticles>(event,"ufsDstar").particles();
      const Particles& particlesDs = apply<UnstableParticles>(event,"ufsDs").particles();
      const Particles& chParticles = apply<ALICE::PrimaryParticles>(event, "app").particles();

      if (CollSystem == "PP5TeV") {

        _c["sow_pp5TeV"]->fill();

        for (const Particle& p : particlesD0) {
          if (p.fromBottom()) continue;
          for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
            _h["den_D0PbPb_pp"+cen]->fill(p.pT()/GeV);
          }
        }

        for (const Particle& p : particlesDplus) {
          if (p.fromBottom()) continue;
          for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
            _h["den_DplusPbPb_pp"+cen]->fill(p.pT()/GeV);
          }
        }

        for (const Particle& p : particlesDstar) {
          if(p.fromBottom()) continue;
          for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
            _h["den_DstarPbPb_pp"+cen]->fill(p.pT()/GeV);
          }
        }

        for (const Particle& p : particlesDs) {
          if (p.fromBottom()) continue;
          for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
            _h["den_DsPbPb_pp"+cen]->fill(p.pT()/GeV);
          }
        }

        for (const Particle& p : chParticles) {
          for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
            _h["den_mult_PbPb_pp"+cen]->fill(p.pT()/GeV);
          }
        }

      }

      if (CollSystem == "PP2TeV") {

        _c["sow_pp2TeV"]->fill();

        for (const Particle& p : particlesD0) {
          if (p.fromBottom()) continue;
          _h["den_D0PbPb_pp2TeV"]->fill(p.pT()/GeV);
        }

        for (const Particle& p : particlesDplus) {
          if (p.fromBottom()) continue;
          _h["den_DplusPbPb_pp2TeV"]->fill(p.pT()/GeV);
        }

        for (const Particle& p : particlesDstar) {
          if (p.fromBottom()) continue;
          _h["den_DstarPbPb_pp2TeV"]->fill(p.pT()/GeV);
        }

      }


      // The centrality projection.
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"V0M");
      const double cent = centProj();
      if(cent >= 80.) vetoEvent;

      if (CollSystem == "PBPB5TeV") {

        string cen("");
        if (cent < 10.)  cen = "00-10";
        else if(cent >= 30. && cent < 50.)  cen = "30-50";
        else if (cent >= 60. && cent < 80.) cen = "60_80";

        if (cen != "") {

          _c["sow_PbPb5TeV"+cen]->fill();

          for (const Particle& p : particlesD0) {
            if (p.fromBottom()) continue;
            _h["D0Pt"+cen]->fill(p.pT()/GeV);
            _h["den_Dplus_D0"+cen]->fill(p.pT()/GeV);
            _h["den_Dstar_D0"+cen]->fill(p.pT()/GeV);
            _h["den_Ds_D0"+cen]->fill(p.pT()/GeV);
            _h["num_D0PbPb_pp"+cen]->fill(p.pT()/GeV);
          }

          for (const Particle& p : particlesDplus) {
            if (p.fromBottom()) continue;
            _h["DplusPt_0010"+cen]->fill(p.pT()/GeV);
            _h["num_Dplus_D0"+cen]->fill(p.pT()/GeV);
            _h["den_Ds_Dplus"+cen]->fill(p.pT()/GeV);
            _h["num_DplusPbPb_pp"+cen]->fill(p.pT()/GeV);
          }

          for (const Particle& p : particlesDstar) {
            if (p.fromBottom()) continue;
            _h["DstarPt_0010"+cen]->fill(p.pT()/GeV);
            _h["num_Dstar_D0"+cen]->fill(p.pT()/GeV);
            _h["num_DstarPbPb_pp"+cen]->fill(p.pT()/GeV);
          }

          for (const Particle& p : particlesDs) {
            if (p.fromBottom())  continue;
            _h["DsPt_0010"+cen]->fill(p.pT()/GeV);
            _h["num_Ds_D0"+cen]->fill(p.pT()/GeV);
            _h["num_Ds_Dplus"+cen]->fill(p.pT()/GeV);
            _h["num_DsPbPb_pp"+cen]->fill(p.pT()/GeV);
          }

          for (const Particle& p : chParticles){
            _h["num_mult_PbPb_pp"+cen]->fill(p.pT()/GeV);
          }
        }

      }
      else if (CollSystem == "PBPB2TeV") {
        if (cent < 10.) {
          _c["sow_PbPb2TeV"]->fill();
          for (const Particle& p : particlesD0) {
            if (p.fromBottom()) continue;
            _h["num_D0PbPb_pp2TeV"]->fill(p.pT()/GeV);
          }
          for (const Particle& p : particlesDplus) {
            if (p.fromBottom()) continue;
            _h["num_DplusPbPb_pp2TeV"]->fill(p.pT()/GeV);
          }
          for (const Particle& p : particlesDstar) {
            if (p.fromBottom()) continue;
            _h["num_DstarPbPb_pp2TeV"]->fill(p.pT()/GeV);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for (auto& item  : _h) {
        for (const string& cen : vector<string>{ "00-10", "30-50", "60-80" }) {
          if (item.first.find("_pp") != string::npos) {
            scale(item.second, _n[cen]/_c["sow_pp5TeV"]->sumW());
            continue;
          }
          if (item.first.find("2TeV") != string::npos) {
            const double sf = (item.first.find("num_") != string::npos)? 1.0 : _n["00-10"];
            scale(item.second, sf / _c["sow_pp2TeV"]->sumW());
            continue;
          }
          if (item.first.find("mult") != string::npos) {
            if (item.first.find("num_") != string::npos) {
              scale(item.second, 1./_c["sow_PbPb5TeV_"+cen]->sumW());
            }
            else {
              scale(item.second, _n[cen] / _c["sow_pp5TeV"]->sumW());
            }
            continue;
          }
          if (item.first.find(cen) != string::npos) {
            const double sf = (item.first.find("PbPb") != string::npos)? 1.0 : 0.5;
            scale(item.second, sf / _c["sow_PbPb5TeV_"+cen]->sumW());
          }
        }
      }

      for (auto& item : _e) {
        if (item.first.find("_") == string::npos)  continue;
        divide(_h["num_"+item.first], _h["den_"+item.first], item.second);
      }

      for (const string& cen : vector<string>{ "00-10", "30-50", "60-80", "2TeV" }) {
        mkAverage("average"+cen, { "D0PbPb_pp"+cen, "DplusPbPb_pp"+cen, "DstarPbPb_pp"+cen });
        if (cen == "2TeV")  continue;
        divide(_e["average"+cen], _e["mult_PbPb_pp"+cen], _e["avgDmult"+cen]);
      }

    }

    ///@}


    /// @name Histograms
    ///@{

    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    map<string, Estimate1DPtr> _e;
    map<string, double> _n{ {"00-10", 1572.}, { "30-50", 264.8 }, { "60-80", 28.31} };

    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2018_I1669819);

}
