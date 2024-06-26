// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> omega phi, omega omega and phi phi
  class BELLE_2012_I1090664 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2012_I1090664);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==223 || Cuts::pid==333), "UFS");
      // counters
      for(unsigned int ix=0; ix<3; ++ix) {
        for(unsigned int iy=0; iy<2; ++iy) {
          book(_sigma[ix][iy],
               "TMP/c_"+toString(ix+1)+"_"+toString(iy+1),
               refData(1, ix+1, 1+iy));
        }
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find any rho mesons
      Particles mesons=apply<UnstableParticles>(event, "UFS").particles();
      for (unsigned int ix=0; ix<mesons.size(); ++ix) {
       	if (mesons[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(mesons[ix],nRes,ncount);
      	bool matched = false;
      	for (unsigned int iy=ix+1; iy<mesons.size(); ++iy) {
      	  if (mesons[iy].children().empty()) continue;
       	  map<long,int> nRes2=nRes;
       	  int ncount2 = ncount;
       	  findChildren(mesons[iy],nRes2,ncount2);
       	  if (ncount2 !=0 ) continue;
       	  matched=true;
      	  for (const auto& val : nRes2) {
      	    if (val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  if (matched) {
            if (mesons[ix].pid()!=mesons[iy].pid()) {
              for (unsigned int iy=0; iy<2; ++iy) {
                _sigma[0][iy]->fill(sqrtS());
              }
            }
            else if (mesons[ix].pid()==333) {
              for (unsigned int iy=0; iy<2; ++iy) {
                _sigma[1][iy]->fill(sqrtS());
              }
            }
            else if (mesons[ix].pid()==223) {
              for (unsigned int iy=0; iy<2; ++iy) {
                _sigma[2][iy]->fill(sqrtS());
              }
            }
            break;
      	  }
      	}
      	if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for (unsigned int ix=0; ix<3; ++ix) {
        for (unsigned int iy=0; iy<2; ++iy) {
          scale(_sigma[ix][iy], fact);
          Estimate1DPtr tmp;
          book(tmp, 1, ix+1, 1+iy);
          barchart(_sigma[ix][iy], tmp);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2012_I1090664);

}
