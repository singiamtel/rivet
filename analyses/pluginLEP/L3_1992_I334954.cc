// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"

namespace Rivet {


  /// @brief Event shapes at MZ
  class L3_1992_I334954 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_1992_I334954);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Sphericity(fs), "Sphericity");
      Thrust thrust(fs);
      declare(thrust, "Thrust"    );
      declare(Hemispheres(thrust), "Hemispheres");
      declare(ParisiTensor(fs), "Parisi");
      declare(FastJets(cfs, JetAlg::DURHAM, 0.7), "DurhamJets");
      declare(FastJets(cfs, JetAlg::JADE  , 0.7), "JadeJets"  );
      // histograms
      book(_histThrust    ,  1, 1, 1);
      book(_histMajor     ,  2, 1, 1);
      book(_histMinor     ,  3, 1, 1);
      book(_histOblateness,  4, 1, 1);
      book(_histJade      ,  6, 1, 1);
      book(_histDurham    ,  7, 1, 1);
      book(_histH3        ,  8, 1, 1);
      book(_histH4        ,  9, 1, 1);
      book(_histSphericity, 10, 1, 1);
      book(_histAplanarity, 11, 1, 1);
      book(_histC         , 12, 1, 1);
      book(_histD         , 13, 1, 1);
      book(_histMJetHeavy , 14, 1, 1);
      book(_histMJetLight , 15, 1, 1);
      book(_histMult      , 16, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // 5 charged particles
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if(cfs.particles().size()<5) vetoEvent;
      // charged particle mult
      const string multi_edge = std::to_string(cfs.particles().size()) + ".0";
      _histMult->fill(multi_edge);
      // Sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _histSphericity->fill(sphericity.sphericity());
      _histAplanarity->fill(sphericity.aplanarity());
      // thrust related
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _histThrust    ->fill(thrust.thrust());
      _histMajor     ->fill(thrust.thrustMajor());
      _histMinor     ->fill(thrust.thrustMinor());
      _histOblateness->fill(thrust.oblateness());
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _histMJetHeavy->fill(hemi.scaledM2high());
      _histMJetLight->fill(hemi.scaledM2low());
      // C and D
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _histC->fill(parisi.C());
      _histD->fill(parisi.D());
      // jet rates
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
      _histDurham->fill(y23);
      const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
      y23 = jadejet.clusterSeq()->exclusive_ymerge_max(2);
      _histJade->fill(y23);
      // fox-wolfram moments
      double H3=0, H4=0;
      const FinalState&  fs = apply<FinalState>(event, "FS");
      for(const Particle & p1 : fs.particles()) {
	double modp1 = p1.p3().mod();
	for(const Particle & p2 : fs.particles()) {
	  double modp2 = p2.p3().mod();
	  double cTheta = p1.p3().dot(p2.p3())/modp1/modp2;
	  double pre = modp1*modp2/sqr(sqrtS());
	  H3 +=   0.5*pre*cTheta*(5.*sqr(cTheta)-3.);
	  H4 += 0.125*pre*((35.*sqr(cTheta)-30.)*sqr(cTheta)+3.);
	}
      }
      _histH3->fill(H3);
      _histH4->fill(H4);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histThrust    , 1./sumOfWeights());
      scale(_histMajor     , 1./sumOfWeights());
      scale(_histMinor     , 1./sumOfWeights());
      scale(_histOblateness, 1./sumOfWeights());
      scale(_histJade      , 1./sumOfWeights());
      scale(_histDurham    , 1./sumOfWeights());
      scale(_histH3        , 1./sumOfWeights());
      scale(_histH4        , 1./sumOfWeights());
      scale(_histSphericity, 1./sumOfWeights());
      scale(_histAplanarity, 1./sumOfWeights());
      scale(_histC         , 1./sumOfWeights());
      scale(_histD         , 1./sumOfWeights());
      scale(_histMJetHeavy , 1./sumOfWeights());
      scale(_histMJetLight , 1./sumOfWeights());
      // percentage and bin width
      scale(_histMult, 100./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _histThrust, _histMajor, _histMinor, _histOblateness;
    Histo1DPtr _histJade,_histDurham;
    Histo1DPtr _histH3,_histH4;
    Histo1DPtr _histSphericity, _histAplanarity;
    Histo1DPtr _histC, _histD;
    Histo1DPtr _histMJetHeavy, _histMJetLight;
    BinnedHistoPtr<string> _histMult;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(L3_1992_I334954);

}
