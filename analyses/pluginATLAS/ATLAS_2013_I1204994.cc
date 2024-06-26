// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {


  /// @brief Upsilon(1,2,3S) production at 7TeV
  class ATLAS_2013_I1204994 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1204994);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(UnstableParticles(Cuts::pid==553 || Cuts::pid==100553 || Cuts::pid==200553), "UFS");

      // branching ratios
      _br = {0.0248,0.0193,0.0218};

      // total cross section
      book(_h_total,1,1,1);

      // double differential
      const vector<double> ybins = {0.,1.2,2.25};
      for(unsigned int ix=0;ix<3;++ix) {
      	book(_h_Upsilon_fid_y[ix],5,1,ix+1);
      	book(_h_Upsilon_y    [ix],9,1,ix+1);
      	book(_h_Upsilon_y_r  [ix],"TMP/Ups"+toString(ix),refData(12,1,1));
        book(_h_Upsilon_fid[ix], ybins);
        book(_h_Upsilon[ix],     ybins);
        book(_h_Upsilon_r[ix],   ybins);
        for (size_t iy=1; iy < ybins.size(); ++iy) {
      	  book(_h_Upsilon_fid[ix]->bin(iy), ix+2,1,iy);
      	  book(_h_Upsilon[ix]->bin(iy), ix+6,1,iy);
      	  book(_h_Upsilon_r[ix]->bin(iy), "TMP/Ups"+toString(ix)+"_"+toString(iy), refData(11,1,iy));
      	}
      }
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           Particles & mup, Particles & mum, unsigned int & ngamma) {
      for (const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::MUON ) {
          ++nstable;
          mum.push_back(p);
        }
        else if (id == PID::ANTIMUON) {
          ++nstable;
          mup.push_back(p);
        }
        else if (id == PID::PI0 || id == PID::K0S || id == PID::K0L ) {
          ++nstable;
        }
        else if (id == PID::GAMMA && p.children().empty() ) {
          ++ngamma;
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, mup, mum,ngamma);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles()) {
      	// check if in fiducal region
      	unsigned int nstable=0, ngamma=0;
      	Particles mup,mum;
      	findDecayProducts(p,nstable,mup,mum,ngamma);
      	bool fiducial = false;
      	if (mup.size()==1 && mum.size()==1 && nstable==ngamma+2) {
      	  fiducial = true;
      	  if(mup[0].perp()<4. || mup[0].abseta()>2.4) fiducial = false;
      	  if(mum[0].perp()<4. || mum[0].abseta()>2.4) fiducial = false;
      	}
      	// pT and rapidity
      	double absrap = p.absrap();
      	double xp = p.perp();
      	// type of upsilon
      	unsigned int iups=p.pid()/100000;
      	if (fiducial) {
      	  if (xp<70.) _h_Upsilon_fid_y[iups]->fill(absrap);
      	  _h_Upsilon_fid[iups]->fill(absrap,xp);
      	}
      	if (xp<70.) {
          _h_Upsilon_y[iups]->fill(absrap);
          _h_Upsilon_y_r[iups]->fill(absrap);
        }
      	_h_Upsilon[iups]->fill(absrap,xp);
      	_h_Upsilon_r[iups]->fill(absrap,xp);
      	if (absrap<2.25 && xp<70) _h_total->fill(_edges[iups],_br[iups]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_total,crossSection()/nanobarn/sumOfWeights());

      for (unsigned int ix=0; ix<3; ++ix) {
      	scale(_h_Upsilon_fid_y[ix], 0.5*crossSection()/picobarn/sumOfWeights());
      	scale(_h_Upsilon_y[ix],     0.5*_br[ix]*crossSection()/picobarn/sumOfWeights());
      	scale(_h_Upsilon_y_r[ix],   0.5*_br[ix]*crossSection()/picobarn/sumOfWeights());

      	scale(_h_Upsilon_fid[ix], 0.5*crossSection()/femtobarn/sumOfWeights());
      	scale(_h_Upsilon[ix],     0.5*_br[ix]*crossSection()/femtobarn/sumOfWeights());
      	scale(_h_Upsilon_r[ix],   0.5*_br[ix]*crossSection()/femtobarn/sumOfWeights());
        divByGroupWidth(_h_Upsilon_fid[ix]);
        divByGroupWidth(_h_Upsilon[ix]);
        divByGroupWidth(_h_Upsilon_r[ix]);
      }

      // ratios
      for (size_t iy=1; iy<_h_Upsilon[0]->numBins()+1; ++iy) {
      	Estimate1DPtr tmp;
      	book(tmp,10,1,iy);
      	divide(_h_Upsilon_r[1]->bin(iy), _h_Upsilon_r[0]->bin(iy), tmp);
      	book(tmp,11,1,iy);
      	divide(_h_Upsilon_r[2]->bin(iy), _h_Upsilon_r[0]->bin(iy), tmp);
      }
      for (unsigned int iy=0; iy<2; ++iy) {
      	Estimate1DPtr tmp;
      	book(tmp,12,1,iy+1);
      	divide(_h_Upsilon_y_r[iy+1], _h_Upsilon_y_r[0], tmp);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_total;
    Histo1DPtr _h_Upsilon_fid_y[3],_h_Upsilon_y[3],_h_Upsilon_y_r[3];
    Histo1DGroupPtr _h_Upsilon[3],_h_Upsilon_fid[3],_h_Upsilon_r[3];
    vector<double> _br;
    vector<string> _edges = { "P P --> UPSILON(1S) < MU+ MU- > X",
                              "P P --> UPSILON(2S) < MU+ MU- > X",
                              "P P --> UPSILON(3S) < MU+ MU- > X" };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1204994);

}
