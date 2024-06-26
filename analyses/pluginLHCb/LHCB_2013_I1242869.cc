// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c at 7 TeV
  class LHCB_2013_I1242869 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1242869);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ichi=0;ichi<3;++ichi) {
	book(_h_chi[ichi],"TMP/h_CHI_"+toString(ichi),refData(1,1,1));
	book(_c_chi[ichi],"TMP/c_CHI_"+toString(ichi));
      }
      book(_r,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_edges.empty()) _edges = _r->xEdges();
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==10441 ||
					     Cuts::pid==20443 ||
					     Cuts::pid==445)) {
	// prompt
	if(p.fromBottom()) continue;
	// J/psi /gamma mode
	if(p.children().size()!=2) continue;
	Particle Jpsi;
	if(p.children()[0].pid()==22 && p.children()[1].pid()==443) {
	  Jpsi=p.children()[1];
	}
	else if(p.children()[1].pid()==22 && p.children()[0].pid()==443) {
	  Jpsi=p.children()[0];
	}
	else
	  continue;
	double absrap=Jpsi.absrap();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int ichi = 0;
	if(p.pid()==20443) ichi=1;
	else if(p.pid()==445) ichi=2;
	double xp=Jpsi.perp();
	_h_chi[ichi]->fill(xp);
	if(xp>4. && xp<20.) {
          _c_chi[ichi]->fill();
          if     (ichi==0) _r->fill(_edges[0]);
          else if(ichi==2) _r->fill(_edges[1]);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // chi_c to gamma J/psi branching ratios from PDG 2021
      vector<double> br = {0.014,0.343,0.190};
      // divide out the branching ratio for mode used to get total rate
      for(unsigned int ichi=0;ichi<3;++ichi) {
	scale(_h_chi[ichi],1./br[ichi]);
	scale(_c_chi[ichi],1./br[ichi]);
      }
      // ratio chi_c2/chi_c2
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      divide(_h_chi[2],_h_chi[1],tmp);
      _r->bin(1).scaleW(double(1./br[0]/ _c_chi[2]->val()));
      _r->bin(2).scaleW(1./br[2]/ _c_chi[1]->val());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_chi[3];
    BinnedHistoPtr<string> _r;
    CounterPtr _c_chi[3];
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1242869);

}
