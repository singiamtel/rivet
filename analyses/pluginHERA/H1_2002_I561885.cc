// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of D*+- meson production in deep inelastic scattering at HERA (H1)
  class H1_2002_I561885 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2002_I561885);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {


      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "Dstars");
      //Cuts::abspid == PID::DSTARPLUS

      // Initialise and register projections


      Histo1DPtr dummy; //Introducing


      // Book histograms
      book(_h["W_GeV"], 2, 1, 1);
      book(_h["p_tD*"], 3, 1, 1);
      book(_h["logx"], 4, 1, 1);
      book(_h["etaD*"], 5, 1, 1);
      book(_h_Q2      , 6, 1, 1);
      book(_h["Z_D*"], 7, 1, 1);

      book(_h_Q2eta, {-1.5, -0.5, 0.5, 1.5}, {"d08-x01-y01", "d08-x01-y02", "d08-x01-y03"});
      book(_h_Q2pt, {1.5, 4., 10.}, {"d09-x01-y01", "d09-x01-y02"});
      
      book(_h_eta1, {0., 0.25, 0.5, 1.}, {"d10-x01-y01", "d10-x01-y02", "d10-x01-y03"});
      book(_h_eta2, {1.5, 2.5, 4., 10.}, {"d11-x01-y01", "d11-x01-y02", "d11-x01-y03"});
      book(_h_zD1, {1.5, 2.5, 4., 10.}, {"d12-x01-y01", "d12-x01-y02", "d12-x01-y03"});
      _axes["Q2"]    = YODA::Axis<double>{1., 2.371, 4.218, 10.0, 17.78, 31.62,100.0};
      _axes["Q2eta"] = YODA::Axis<double>{1.0, 9.0, 20.0,100.0};
      _axes["Q2pt"]  = YODA::Axis<double>{1.0, 3.0, 8.0, 20.0, 100.0};
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty()) {
        _edges["Q2"   ] = _h_Q2->xEdges();
        _edges["Q2eta"] = _h_Q2eta->bins()[1]->xEdges();
        _edges["Q2pt" ] = _h_Q2pt ->bins()[1]->xEdges();
      }
      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");

       // Q2 and inelasticity cuts
      if (!inRange(kin.Q2(), 1.0*GeV2, 100*GeV2)) vetoEvent;
      if (!inRange(kin.y(), 0.05, 0.7)) vetoEvent;

        // D* reconstruction
      // const Particles unstables = apply<ParticleFinder>(event, "Dstars").particles(Cuts::pT > 1.5*GeV && Cuts::abseta < 1.5);
      const Particles unstables = apply<ParticleFinder>(event, "Dstars").particles(Cuts::pT > 1.5*GeV );
      const Particles dstars = select(unstables, [](const Particle& p){ return p.abspid() == PID::DSTARPLUS; });
      if (dstars.empty()) vetoEvent;
      // MSG_DEBUG("#D* = " << dstars.size());
      //const Particle& dstar = dstars.front();
      for (const Particle& dstar : dstars){
        const double zD = (dstar.E() - dstar.pz()) / (2*kin.beamLepton().E()*kin.y());


        // Single-differential histograms
        _h["p_tD*"]->fill(dstar.pT()/GeV);
        _h["etaD*"]->fill(dstar.eta());
        _h["Z_D*"]->fill(zD/GeV);
        
        size_t idx = _axes["Q2"].index(kin.Q2()/GeV2);
        string edge = "OTHER";
        if(idx && idx <= _edges["Q2"].size()) edge=_edges["Q2"][idx-1];
        _h_Q2->fill(edge);
        
        _h["W_GeV"]->fill(sqrt(kin.W2())/GeV);
        _h["logx"]->fill(log10(kin.x()));

        // Double-differential (y,Q2) histograms

        idx = _axes["Q2eta"].index(kin.Q2()/GeV2);
        edge = "OTHER";
        if(idx && idx <= _edges["Q2eta"].size()) edge=_edges["Q2eta"][idx-1];
        _h_Q2eta->fill(dstar.eta(),edge);
        idx = _axes["Q2pt"].index(kin.Q2()/GeV2);
        edge = "OTHER";
        if(idx && idx <= _edges["Q2pt"].size()) edge=_edges["Q2pt"][idx-1];
        _h_Q2pt->fill(dstar.pT(),edge);
        _h_eta1->fill(zD, dstar.eta());
        _h_eta2->fill(dstar.pT(), dstar.eta());
        _h_zD1->fill(dstar.pT(), zD/GeV);

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


     const double sf = crossSection()/nanobarn/sumOfWeights();
      scale(_h["p_tD*"], sf);
      scale(_h["etaD*"], sf);
      scale(_h["Z_D*"], sf);
      scale(_h_Q2     , sf);
      for(auto & b : _h_Q2->bins()) {
        const size_t idx=b.index();
        b.scaleW(1./_axes["Q2"].width(idx));
      }
      scale(_h["W_GeV"], sf);
      scale(_h["logx"], sf);

      scale(_h_Q2eta, sf);
      for(auto & histo : _h_Q2eta->bins()) {
        for(auto & b : histo->bins()) {
          const size_t idx=b.index();
          b.scaleW(1./_axes["Q2eta"].width(idx));
        }
      }
      scale(_h_Q2pt, sf);
      for(auto & histo : _h_Q2pt->bins()) {
        for(auto & b : histo->bins()) {
          const size_t idx=b.index();
          b.scaleW(1./_axes["Q2pt"].width(idx));
        }
      }
      scale(_h_eta1, sf);
      scale(_h_eta2, sf);
      scale(_h_zD1, sf);
      // distributionsa are really double differential divide by group width
      _h_Q2eta->divByGroupWidth();
      _h_Q2pt ->divByGroupWidth();
      _h_eta1 ->divByGroupWidth();
      _h_eta2 ->divByGroupWidth();
      _h_zD1  ->divByGroupWidth();
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    HistoGroupPtr<double,string> _h_Q2eta, _h_Q2pt;
    Histo1DGroupPtr _h_eta1,_h_eta2,_h_zD1;
    BinnedHistoPtr<string> _h_Q2;
    map<string,YODA::Axis<double> > _axes;
    map<string,vector<string> > _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_2002_I561885);

}
