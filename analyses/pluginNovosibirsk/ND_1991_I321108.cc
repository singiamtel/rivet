// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ND measurement of exclusive hadronic final states
  class ND_1991_I321108 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ND_1991_I321108);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book hists
      for(unsigned int ix=0;ix<4;++ix)
        book(_nOmegaPi[ix], 1+ix,1,1);
      book(_nOmegaPi[4], 10,1,3);
      book(_n2Pi    ,11,1,1);
      book(_n3Pi[0] , 5,1,1);
      book(_n3Pi[1] ,10,1,4);
      book(_n4PiC[0], 7,1,1);
      book(_n4PiC[1],10,1,1);
      book(_n4PiN[0], 8,1,1);
      book(_n4PiN[1],10,1,2);
      book(_nEtaPiPi, 6,1,1);
      book(_nKC     ,12,1,1);
      book(_nKN     ,13,1,1);
      book(_n5Pi    ,14,1,1);
      vector<string> energies({"0.661", "0.671", "0.681", "0.691", "0.701", "0.711", "0.725", "0.735", "0.745", "0.755",
          "0.765", "0.805", "0.81", "0.815", "0.825", "0.83", "0.831", "0.841", "0.85", "0.851", "0.861", "0.87", "0.871",
          "0.881", "0.89", "0.891", "0.901", "0.91", "0.911", "0.921", "0.93", "0.943", "0.95", "0.955", "0.956", "0.963",
          "0.97", "0.973", "0.98", "0.983", "0.99", "0.993", "1.0", "1.003", "1.005", "1.02", "1.036", "1.04", "1.05",
          "1.059", "1.06", "1.07", "1.075", "1.08", "1.09", "1.099", "1.1", "1.11", "1.12", "1.13", "1.139", "1.14", "1.15",
          "1.16", "1.17", "1.179", "1.18", "1.19", "1.2", "1.21", "1.219", "1.22", "1.23", "1.24", "1.25", "1.259", "1.26",
          "1.27", "1.28", "1.29", "1.299", "1.3", "1.31", "1.32", "1.325", "1.33", "1.339", "1.34", "1.35", "1.36", "1.37",
          "1.375", "1.379", "1.38", "1.39", "1.395", "1.4", "1.5", "1.6"});
      for(const string& en : energies) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          ecms = en;
          break;
        }
      }
      if(ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==2) {
	if(nCount[211]==1&&nCount[-211]==1)
	  _n2Pi->fill(ecms);
	if(nCount[321]==1&&nCount[-321]==1)
	  _nKC->fill(ecms);
	if(nCount[310]==1&&nCount[130]==1)
	  _nKN->fill(ecms);
      }
      else if(ntotal==3) {
	if(nCount[211]==1&&nCount[-211]==1&&nCount[111]==1) {
	  _n3Pi[0]->fill(ecms);
	  _n3Pi[1]->fill(ecms);
        }
      }
      else if(ntotal==4) {
	if(nCount[211]==2&&nCount[-211]==2) {
	  _n4PiC[0]->fill(ecms);
	  _n4PiC[1]->fill(ecms);
        }
	else if(nCount[211]==1&&nCount[-211]==1&&nCount[111]==2) {
	  _n4PiN[0]->fill(ecms);
	  _n4PiN[1]->fill(ecms);
        }
      }
      else if(ntotal==5) {
	if(nCount[211]==2&&nCount[-211]==2&&nCount[111]==1)
	  _n5Pi->fill(ecms);
      }

      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the eta
	if(p.pid()==221) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
	  if(ncount!=2) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211) {
	      if(val.second !=1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched)
	    _nEtaPiPi->fill(ecms);
	}
	else if(p.pid()==223) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
	  if(ncount!=1) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==111) {
	      if(val.second !=1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
            for(unsigned int ix=0;ix<5;++ix)
              _nOmegaPi[ix]->fill(ecms);
          }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
        scale(_nOmegaPi[ix], crossSection()/ sumOfWeights() /nanobarn);
      scale(_n2Pi    , crossSection()/ sumOfWeights() /nanobarn);
      scale(_nEtaPiPi, crossSection()/ sumOfWeights() /nanobarn);
      scale(_nKC     , crossSection()/ sumOfWeights() /nanobarn);
      scale(_nKN     , crossSection()/ sumOfWeights() /nanobarn);
      scale(_n5Pi    , crossSection()/ sumOfWeights() /nanobarn);
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_n3Pi[ix] , crossSection()/ sumOfWeights() /nanobarn);
        scale(_n4PiC[ix], crossSection()/ sumOfWeights() /nanobarn);
        scale(_n4PiN[ix], crossSection()/ sumOfWeights() /nanobarn);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nOmegaPi[5],_n2Pi,_n3Pi[2],_n4PiC[2],_n4PiN[2],_nEtaPiPi,_nKC,_nKN,_n5Pi;
    string ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ND_1991_I321108);


}
