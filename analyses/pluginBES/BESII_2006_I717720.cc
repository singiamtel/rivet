// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief e+ e- > hadrons
  class BESII_2006_I717720 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2006_I717720);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==411 ||
				Cuts::abspid==421 ), "UFS");
      declare(FinalState(), "FS");

      // Book histograms
      book(_sigma[0], 1,1,1);
      for(unsigned int ix=1;ix<4;++ix)
	book(_sigma[ix], 2,1,ix);
      string eCent = getOption<string>("ECENT", std::to_string(sqrtS()/GeV));
      double ee = std::stod(eCent);
      for(unsigned int ix=0;ix<3;++ix) {
        for(const auto & en : _sigma[ix].binning().edges<0>()) {
          if(fuzzyEquals(std::stod(en),ee,1e-5)) {
            _ecms[ix]=en;
            break;
          }
        }
      }
      if(_ecms[0].empty() && _ecms[1].empty() && _ecms[2].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
      _ecms[3]=_ecms[2];
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
      	vetoEvent;
      }
      else if(nCount[-11]==1 and nCount[11]==1 && ntotal==2+nCount[22]) {
      	vetoEvent;
      }
      // everything else
      else {
      	if(!_ecms[0].empty()) _sigma[0]->fill(_ecms[0]);
	if(!ufs.particles().empty()) {
	  if(!_ecms[1].empty()) _sigma[1]->fill(_ecms[1]);
	  if(ufs.particles()[0].abspid()==421) {
	    if(!_ecms[2].empty()) _sigma[2]->fill(_ecms[2]);
          }
	  else {
	    if(!_ecms[3].empty()) _sigma[3]->fill(_ecms[3]);
          }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/sumOfWeights()/nanobarn;
      for(unsigned int ix=0;ix<4;++ix)
        scale(_sigma[ix],fact);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[4];
    string _ecms[4];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESII_2006_I717720);

}
