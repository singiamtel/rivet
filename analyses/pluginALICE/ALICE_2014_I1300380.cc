// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


  class ALICE_2014_I1300380 : public Analysis {
  public:

    ALICE_2014_I1300380()
      : Analysis("ALICE_2014_I1300380")
    {}


  public:

    void init() {
      const UnstableParticles cfs(Cuts::absrap<0.5);
      declare(cfs, "CFS");

      // Plots from the paper
      book(_histPtSigmaStarPlus        ,"d01-x01-y01");    // Sigma*+
      book(_histPtSigmaStarMinus       ,"d01-x01-y02");    // Sigma*-
      book(_histPtSigmaStarPlusAnti    ,"d01-x01-y03");    // anti Sigma*-
      book(_histPtSigmaStarMinusAnti   ,"d01-x01-y04");    // anti Sigma*+
      book(_histPtXiStar               ,"d02-x01-y01");    // 0.5 * (xi star + anti xi star)
      book(_histAveragePt              ,"d03-x01-y01");  // <pT> profile
    }


    void analyze(const Event& event) {
      if(_edges.empty()) _edges=_histAveragePt->xEdges();
      const UnstableParticles& cfs = apply<UnstableParticles>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
	// protections against mc generators decaying long-lived particles
	if ( !(p.hasAncestorWith(Cuts::pid == 310)  || p.hasAncestorWith(Cuts::pid == -310)   || // K0s
	       p.hasAncestorWith(Cuts::pid == 130)  || p.hasAncestorWith(Cuts::pid == -130)   ||     // K0l
	       p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322)  ||     // Xi0
	       p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122)  ||     // Lambda
	       p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222)  ||     // Sigma+/-
	       p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312)  ||     // Xi-/+
	       p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334)  ))     // Omega-/+
	{
	  int aid = p.abspid();
	  if      (aid == 211 ) _histAveragePt->fill(_edges[0], p.pT()/GeV); // pi+
          else if (aid == 321 ) _histAveragePt->fill(_edges[1], p.pT()/GeV); // K+
          else if (aid == 313 ) _histAveragePt->fill(_edges[2], p.pT()/GeV);   // K*(892)0
          else if (aid == 2212) _histAveragePt->fill(_edges[3], p.pT()/GeV);    // proton
          else if (aid == 333 ) _histAveragePt->fill(_edges[4], p.pT()/GeV);    // phi(1020)
	} // end if "rejection of long-lived particles"
        switch (p.pid()) {
	  case 3224:
	    _histPtSigmaStarPlus->fill(p.pT()/GeV);
	    _histAveragePt->fill(_edges[6], p.pT()/GeV);
	    break;
	  case -3224:
	    _histPtSigmaStarPlusAnti->fill(p.pT()/GeV);
	    _histAveragePt->fill(_edges[6], p.pT()/GeV);
	    break;
	  case 3114:
	    _histPtSigmaStarMinus->fill(p.pT()/GeV);
	    _histAveragePt->fill(_edges[6], p.pT()/GeV);
	    break;
	  case -3114:
	    _histPtSigmaStarMinusAnti->fill(p.pT()/GeV);
	    _histAveragePt->fill(_edges[6], p.pT()/GeV);
	    break;
	  case 3324:
	    _histPtXiStar->fill(p.pT()/GeV);
	    _histAveragePt->fill(_edges[7], p.pT()/GeV);
	    break;
	  case -3324:
	    _histPtXiStar->fill(p.pT()/GeV);
	    _histAveragePt->fill(_edges[7], p.pT()/GeV);
	    break;
	  case 3312:
	    _histAveragePt->fill(_edges[5], p.pT()/GeV);
	    break;
	  case -3312:
	    _histAveragePt->fill(_edges[5], p.pT()/GeV);
	    break;
	  case 3334:
	    _histAveragePt->fill(_edges[8], p.pT()/GeV);
	    break;
	  case -3334:
	    _histAveragePt->fill(_edges[8], p.pT()/GeV);
	    break;
        }
      }
    }


    void finalize() {
      scale(_histPtSigmaStarPlus,       1./sumOfWeights());
      scale(_histPtSigmaStarPlusAnti,   1./sumOfWeights());
      scale(_histPtSigmaStarMinus,      1./sumOfWeights());
      scale(_histPtSigmaStarMinusAnti,  1./sumOfWeights());
      scale(_histPtXiStar,              1./sumOfWeights()/ 2.);
    }


  private:
    // plots from the paper
    Histo1DPtr   _histPtSigmaStarPlus;
    Histo1DPtr   _histPtSigmaStarPlusAnti;
    Histo1DPtr   _histPtSigmaStarMinus;
    Histo1DPtr   _histPtSigmaStarMinusAnti;
    Histo1DPtr   _histPtXiStar;
    BinnedProfilePtr<string> _histAveragePt;
    vector<string> _edges;
  };


  RIVET_DECLARE_PLUGIN(ALICE_2014_I1300380);

}
