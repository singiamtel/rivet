// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> J/psi/psi(2S) K+ pi+ pi-
  class BELLE_2010_I871475 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2010_I871475);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(100443);
      BP.addStable(443);
      declare(BP, "BP");
      // histograms
      vector<double> bins[2] = {{0.60,1.46,2.32,3.18,4.04,4.90},
                                {0.60,1.00,1.40,1.80,2.20,2.60}};
      for (size_t ix=0; ix<2; ++ix) {
        book(_c[ix],"TMP/c_"+toString(ix+1));
        for (size_t iy=0; iy<3; ++iy) {
          book(_h_all[ix][iy], 1+ix, 1, 1+iy);
          if (iy>1) continue;
          book(_h_slice[ix][iy], bins[ix]);
          for (auto& b : _h_slice[ix][iy]->bins()) {
            book(b, 3+ix, 1+iy, b.index());
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1}, { 211,1}, {-211,1}, {    443,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1}, { 211,1}, {-211,1}, {    443,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 321,1}, { 211,1}, {-211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-321,1}, { 211,1}, {-211,1}, { 100443,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0;ix<BP.decaying().size();++ix) {
        int imode=-1, sign = 1;
        if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,4,mode1)) {
          sign=1;
          imode=0;
        }
        else if  (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,4,mode1CC)) {
          sign=-1;
          imode=0;
        }
        else if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,4,mode2)) {
          sign=1;
          imode=1;
        }
        else if  (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,4,mode2CC)) {
          sign=-1;
          imode=1;
        }
        else {
          continue;
        }
      	const Particle& Kp  = BP.decayProducts()[ix].at( 321*sign)[0];
       	const Particle& pip = BP.decayProducts()[ix].at( 211*sign)[0];
       	const Particle& pim = BP.decayProducts()[ix].at(-211*sign)[0];
        double m2Kpipi = (Kp .momentum()+pip.momentum()+pim.momentum()).mass2();
        double m2Kpi   = (Kp .momentum()+pim.momentum()).mass2();
        double m2pipi  = (pip.momentum()+pim.momentum()).mass2();
        _h_all[imode][0]->fill(m2Kpipi);
        _h_all[imode][1]->fill(m2Kpi  );
        _h_all[imode][2]->fill(m2pipi );
        _h_slice[imode][0]->fill(m2Kpipi, m2Kpi );
        _h_slice[imode][1]->fill(m2Kpipi, m2pipi);
        _c[imode]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t ix=0; ix<2; ++ix) {
        for (size_t iy=0; iy<3; ++iy) {
          scale(_h_all[ix][iy], 1./ *_c[ix]);
          if(iy>1) continue;
          scale(_h_slice[ix][iy], 1./ *_c[ix]);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_all[2][3];
    Histo1DGroupPtr _h_slice[2][2];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2010_I871475);

}
