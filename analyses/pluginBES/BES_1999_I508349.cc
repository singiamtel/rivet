// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief charm cross sections
  class BES_1999_I508349 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1999_I508349);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(),"UFS");
      for(unsigned int ix=0;ix<4;++ix)
        book(_nD[ix], 1, 1, 1+ix);
      
      if (isCompatibleWithSqrtS(4.03*GeV,1e-3)) {
        book(_h_D0,2,1,1);
        book(_h_Dp,2,1,2);
        _ecms = "4.03";
      }
      else if (isCompatibleWithSqrtS(4.14*GeV,1e-3)) {
        book(_h_D0,3,1,1);
        book(_h_Dp,3,1,2);
        _ecms = "4.14";
      } else {
        MSG_ERROR("Beam energy not supported!");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double nD0(0),nDp(0),nDs(0);
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==411 or Cuts::abspid==421 or Cuts::abspid==431)) {
        if (p.abspid()==421) {
          _h_D0->fill(p.momentum().p3().mod());
          ++nD0;
        }
        else if (p.abspid()==411) {
          _h_Dp->fill(p.momentum().p3().mod());
          ++nDp;
        }
        else {
          ++nDs;
        }
      }
      _nD[0]->fill(_ecms,0.5*nD0);
      _nD[1]->fill(_ecms,0.5*nDp);
      _nD[2]->fill(_ecms,0.5*nDs);
      _nD[3]->fill(_ecms,0.5*(nD0+nDp+nDs));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/sumOfWeights()/nanobarn;
      scale(_h_D0,0.5*fact);
      scale(_h_Dp,0.5*fact);
      for(unsigned int ix=0;ix<4;++ix)
        scale(_nD[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nD[4];
    Histo1DPtr _h_D0,_h_Dp;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BES_1999_I508349);

}
