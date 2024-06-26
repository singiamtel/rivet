// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s -> K+K-pi+ /K+K+pi-
  class BABAR_2010_I878120 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I878120);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      declare(DS, "DS");
      // histos
      book(_h_KK[0],1,1,1);
      book(_h_KK[1],1,1,2);
      book(_h_Kmpi ,1,1,3);
      book(_h_Kppi ,1,1,4);
      book(_h_Kppi2,1,1,5);
      book(_dalitz, "dalitz",50,0.3,3.5,50,0.07,2.5);
      book(_dalitz2, "dalitz2",50,0.3,2.3,50,0.03,2.3);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 211,1}, { 321,1}, {-321,1} };
      static const map<PdgId,unsigned int> & mode1CC = { {-211,1}, { 321,1}, {-321,1} };
      static const map<PdgId,unsigned int> & mode2   = { {-211,1}, { 321,2} };
      static const map<PdgId,unsigned int> & mode2CC = { { 211,1}, {-321,2} };
      // Loop over D+ mesons
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	int sign = DS.decaying()[ix].pid()/DS.decaying()[ix].abspid();
       	if     ( DS.modeMatches(ix,3,mode1  ) ||
		 DS.modeMatches(ix,3,mode1CC)) {
          const Particle & pip = DS.decayProducts()[ix].at( sign*211)[0];
          const Particle & Km  = DS.decayProducts()[ix].at(-sign*321)[0];
          const Particle & Kp  = DS.decayProducts()[ix].at( sign*321)[0];
          double mplus  = (Kp.momentum()+pip.momentum()).mass2();
          double mminus = (Km.momentum()+pip.momentum()).mass2();
          double mKK    = (Kp.momentum()+Km .momentum()).mass2();
          _h_KK[0]->fill(mKK);
          _h_KK[1]->fill(mKK);
          _h_Kppi->fill(mplus);
          _h_Kmpi->fill(mminus);
          _dalitz->fill(mKK,mminus);
	}
	else if( DS.modeMatches(ix,3,mode2  ) ||
		 DS.modeMatches(ix,3,mode2CC)) {
	  const Particle  & pim = DS.decayProducts()[ix].at(-sign*211)[0];
	  const Particles & Kp  = DS.decayProducts()[ix].at( sign*321);
	  double mKpi[2];
	  for(unsigned int ix=0;ix<2;++ix){
	    mKpi[ix] = (Kp[ix].momentum()+pim.momentum()).mass2();
	    _h_Kppi2->fill(sqrt(mKpi[ix]));
	  }
	  _dalitz2->fill(mKpi[0],mKpi[1]);
	  _dalitz2->fill(mKpi[1],mKpi[0]);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK[0]);
      normalize(_h_KK[1],1.,false);
      normalize(_h_Kmpi);
      normalize(_h_Kppi);
      normalize(_h_Kppi2);
      normalize(_dalitz);
      normalize(_dalitz2);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK[2],_h_Kmpi,_h_Kppi,_h_Kppi2;
    Histo2DPtr _dalitz,_dalitz2;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2010_I878120);

}
