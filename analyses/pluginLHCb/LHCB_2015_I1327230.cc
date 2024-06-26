// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B_c production
  class LHCB_2015_I1327230 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2015_I1327230);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ib=0;ib<2;++ib) {
	book(_h_pT[ib], "TMP/h_pT_"+toString(ib), refData(2,1,1));
	book(_h_y [ib], "TMP/h_y_" +toString(ib), refData(3,1,1));
        book(_h_B[ib],{2.0,2.9,3.3,4.5});
	for(unsigned int iy=0;iy<3;++iy) {
          book(_h_B[ib]->bin(iy+1),"TMP/hB_"+toString(ib)+"_"+toString(iy),refData(1,1,1+iy));
	}
	book(_c_Bc[ib],"TMP/c_Bc_"+toString(ib));
      }
      _brBp = getOption<double>("BPBR",1.02e-3);
      _brBc = getOption<double>("BCBR",-1.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==521 or Cuts::pid==541 )) {
	if(p.pid()==541) {
	  if(p.children().size()==2 && 
	     ( (p.children()[0].pid()==211 && p.children()[1].pid()==443) ||
	       (p.children()[1].pid()==211 && p.children()[0].pid()==443) ))
	    _c_Bc[0]->fill(1.);
	  _c_Bc[1]->fill(1.);
	}
        double absrap = p.absrap();
	if(absrap<2. || absrap>4.5) continue;
        double pT = p.perp();
	// select B0
	if(p.pid()==521) {
	  _h_B [1]->fill(absrap,pT);
	  _h_pT[1]->fill(pT);
	  if(pT<20.) _h_y [1]->fill(absrap);
	}
	// select B_c+
	else {
	  _h_B [0]->fill(absrap,pT);
	  _h_pT[0]->fill(pT);
	  if(pT<20.) _h_y [0]->fill(absrap);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // branching ratio Bc -> J/psi pi+
      if (_brBc<0.) {
        Estimate0DPtr tmp2;
        book(tmp2,"TMP/br");
        divide(_c_Bc[0],_c_Bc[1],tmp2);
        _brBc = tmp2->val();
      }
      // scale by B+ br (for B_c directly selected before filling histos)
      scale(_h_B[0], _brBc);
      scale(_h_B[1],_brBp);
      for(unsigned int ix=0;ix<2;++ix) divByGroupWidth(_h_B[ix]);
      scale(_h_pT[0], _brBc);
      scale(_h_pT[1],_brBp);
      scale(_h_y [0], _brBc);
      scale(_h_y [1], _brBp);
      for(unsigned int iy=0;iy<3;++iy) {
	Estimate1DPtr tmp;
	book(tmp,1,1,1+iy);
	divide(_h_B[0]->bin(iy+1),_h_B[1]->bin(iy+1),tmp);
	// convert ratio to %
	tmp->scale(100.);
      }
      Estimate1DPtr tmp;
      book(tmp,2,1,1);
      divide(_h_pT[0],_h_pT[1],tmp);
      // convert ratio to %
      tmp->scale(100.);
      book(tmp,3,1,1);
      divide(_h_y [0],_h_y [1],tmp);
      // convert ratio to %
      tmp->scale(100.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_B[2];
    Histo1DPtr _h_pT[2],_h_y[2];
    CounterPtr _c_Bc[2];
    double _brBp,_brBc;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2015_I1327230);

}
