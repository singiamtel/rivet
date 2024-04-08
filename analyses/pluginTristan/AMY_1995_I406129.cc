// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JadePlugin.hh"

namespace fastjet {

class P_scheme : public JetDefinition::Recombiner {
 public:
  std::string description() const {return "";}
  void recombine(const PseudoJet & pa, const PseudoJet & pb,
		 PseudoJet & pab) const {
    PseudoJet tmp = pa + pb;
    double E = sqrt(tmp.px()*tmp.px() + tmp.py()*tmp.py() + tmp.pz()*tmp.pz());
    pab.reset_momentum(tmp.px(), tmp.py(), tmp.pz(), E);
  }
  void preprocess(PseudoJet & p) const {
    double E = sqrt(p.px()*p.px() + p.py()*p.py() + p.pz()*p.pz());
    p.reset_momentum(p.px(), p.py(), p.pz(), E);
  }
  ~P_scheme() { }
};

class E0_scheme : public JetDefinition::Recombiner {
 public:
  std::string description() const {return "";}
  void recombine(const PseudoJet & pa, const PseudoJet & pb,
		 PseudoJet & pab) const {
    PseudoJet tmp = pa + pb;
    double fact = tmp.E()/sqrt(tmp.px()*tmp.px()+tmp.py()*tmp.py()+tmp.pz()*tmp.pz());
    pab.reset_momentum(fact*tmp.px(), fact*tmp.py(), fact*tmp.pz(), tmp.E());
  }
  void preprocess(PseudoJet & p) const {
    double fact = p.E()/sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz());

    p.reset_momentum(fact*p.px(), fact*p.py(), fact*p.pz(), p.E());
  }
  ~E0_scheme() { }
};

}


namespace Rivet {


  /// @brief AMY jets at
  class AMY_1995_I406129 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(AMY_1995_I406129);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      FinalState fs;
      declare(fs, "FS");
      // Book histograms
      book(_h["jade_P"],  2, 1, 1);
      book(_h["jade_E"],  3, 1, 1);
      book(_h["jade_E0"], 6, 1, 1);
      book(_h["durham"],  4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles particles =  apply<FinalState>(event, "FS").particles();
      MSG_DEBUG("Num particles = " << particles.size());
      PseudoJets pjs;
      double mpi=.13957;
      for (const Particle & p : particles) {
        Vector3 mom = p.p3();
        double energy = p.E();
        if(PID::isCharged(p.pid())) {
          energy = sqrt(mom.mod2()+sqr(mpi));
        }
        else {
          double fact = energy/mom.mod();
          mom *=fact;
        }
        pjs.push_back(fastjet::PseudoJet(mom.x(),mom.y(),mom.z(),energy));
      }
      // durham
      fastjet::JetDefinition durDef(fastjet::ee_kt_algorithm, fastjet::E_scheme);
      fastjet::ClusterSequence durham(pjs,durDef);
      double y_23 = durham.exclusive_ymerge_max(2);
      size_t idx = durhamAxis.index(y_23);
      string label = idx>0 && idx <= _h["durham"]->xEdges().size() ?  _h["durham"]->xEdges()[idx-1] : "OTHER";
      _h["durham"]->fill(label);
      // jade e-scheme
      fastjet::JetDefinition::Plugin *plugin = new fastjet::JadePlugin();
      fastjet::JetDefinition jadeEDef(plugin);
      jadeEDef.set_recombination_scheme(fastjet::E_scheme);
      fastjet::ClusterSequence jadeE(pjs,jadeEDef);
      y_23 = jadeE.exclusive_ymerge_max(2);
      idx = jadeAxis.index(y_23);
      label = idx>0 && idx <= _h["jade_E"]->xEdges().size() ?  _h["jade_E"]->xEdges()[idx-1] : "OTHER";
      _h["jade_E"]->fill(label);
      // jade p-scheme
      fastjet::P_scheme p_scheme;
      fastjet::JetDefinition jadePDef(plugin);
      jadePDef.set_recombiner(&p_scheme);
      fastjet::ClusterSequence jadeP(pjs,jadePDef);
      y_23 = jadeP.exclusive_ymerge_max(2);
      idx = jadeAxis.index(y_23);
      label = idx>0 && idx <= _h["jade_P"]->xEdges().size() ?  _h["jade_P"]->xEdges()[idx-1] : "OTHER";
      _h["jade_P"]->fill(label);
      // jade E0-scheme
      fastjet::E0_scheme e0_scheme;
      fastjet::JetDefinition jadeE0Def(plugin);
      jadeE0Def.set_recombiner(&e0_scheme);
      fastjet::ClusterSequence jadeE0(pjs,jadeE0Def);
      y_23 = jadeE0.exclusive_ymerge_max(2);
      idx = jadeAxis.index(y_23);
      label = idx>0 && idx <= _h["jade_E0"]->xEdges().size() ?  _h["jade_E0"]->xEdges()[idx-1] : "OTHER";
      _h["jade_E0"]->fill(label);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1./sumOfWeights());
      for( auto & hist : _h) {
        for(auto & b: hist.second->bins()) {
          const size_t idx = b.index();
          b.scaleW(hist.first=="durham" ? 1./durhamAxis.width(idx) : 1./jadeAxis.width(idx));
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,BinnedHistoPtr<string>> _h;
    YODA::Axis<double> jadeAxis{0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.0825, 0.1,
                                0.12, 0.14, 0.16, 0.1825, 0.21, 0.24, 0.27, 0.3};
    YODA::Axis<double> durhamAxis{0.0, 0.002, 0.004, 0.0065, 0.0115, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                                  0.0825, 0.1, 0.12, 0.14, 0.16, 0.1825, 0.21, 0.24, 0.27, 0.3};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(AMY_1995_I406129);


}
