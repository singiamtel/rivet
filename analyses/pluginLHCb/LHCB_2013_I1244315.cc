// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi polarization at 7 TeV
  class LHCB_2013_I1244315 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1244315);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // histograms
      // cross section
      _ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      book(_h_Jpsi,_ybins);
      for(unsigned int iy=0;iy<5;++iy) {
	book(_h_Jpsi->bin(iy+1),4,1,1+iy);
	for(unsigned int ix=0;ix<2;++ix) {
	  for(unsigned int iz=0;iz<3;++iz) {
	    string name="TMP/POL_"+toString(ix)+"_"+toString(iy)+"_"+toString(iz);
	    book(_p_Jpsi[ix][iy][iz],name,refData(1+ix,1,1+iy));
	  }
	}
      }
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,  
                           Particles & mup, Particles & mum) {
      for(const Particle & p : mother.children()) {
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
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, mup, mum);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the beams
      const ParticlePair & beams = apply<Beam>(event, "Beams").beams();
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
	// prompt
	if(p.fromBottom()) continue;
      	// pT and rapidity
      	double rapidity = p.rapidity();
      	double xp = p.perp();
	// cross section
	_h_Jpsi->fill(rapidity,xp);
	// polarization
      	unsigned int nstable=0;
      	Particles mup,mum;
      	findDecayProducts(p,nstable,mup,mum);
      	if(mup.size()!=1 || mum.size()!=1 || nstable!=2) continue;
	// find the rapidity interval
	if(rapidity<=2. || rapidity>=4.5) continue;
	unsigned int iy=0;
	for(iy=0;iy<5;++iy) if(rapidity<_ybins[iy+1]) break;
	// first the CS frame
	// first boost so upslion momentum =0 in z direction
	Vector3 beta = p.momentum().betaVec();
	beta.setX(0.);beta.setY(0.);
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(beta);
	FourMomentum pp = boost.transform(p.momentum());
	// and then transverse so pT=0
        beta = pp.betaVec();
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(beta);
	// get all the momenta in this frame
	Vector3 muDirn = boost2.transform(boost.transform(mup[0].momentum())).p3().unit();
	FourMomentum p1 = boost2.transform(boost.transform(beams. first.momentum()));
	FourMomentum p2 = boost2.transform(boost.transform(beams.second.momentum()));
	if(beams.first.momentum().z()<0.) swap(p1,p2);
	if(p.rapidity()<0.) swap(p1,p2);
	Vector3 axisy = (p1.p3().cross(p2.p3())).unit();
	Vector3 axisz(0.,0.,1.);
	Vector3 axisx = axisy.cross(axisz);
	double cTheta = axisz.dot(muDirn);
	double cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_Jpsi[1][iy][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Jpsi[1][iy][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Jpsi[1][iy][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	// now for the HX frame
	beta = p.momentum().betaVec();
	boost = LorentzTransform::mkFrameTransformFromBeta(beta);
	axisz = pp.p3().unit();
	axisx = axisy.cross(axisz);
	cTheta = axisz.dot(muDirn);
	cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_Jpsi[0][iy][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Jpsi[0][iy][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Jpsi[0][iy][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/nanobarn/sumOfWeights();
      scale(_h_Jpsi,factor);
      divByGroupWidth(_h_Jpsi);
      for(unsigned int iy=0;iy<5;++iy) {
	// Loop over frame definition
	for(unsigned int iframe=0;iframe<2;++iframe) {
          // book scatters
          Estimate1DPtr lTheta,lPhi,lThetaPhi;
          book(lTheta   ,iframe+1,1,1+iy);
          book(lPhi     ,iframe+1,3,1+iy);
          book(lThetaPhi,iframe+1,2,1+iy);
	  // histos for the moments
	  Profile1DPtr moment[3];
	  for(unsigned int ix=0;ix<3;++ix)
	    moment[ix] = _p_Jpsi[iframe][iy][ix];
	  // loop over bins
	  for(unsigned int ibin=1;ibin<=moment[0]->bins().size();++ibin) {
            // extract moments and errors
            double val[3],err[3];
            // m1 = lTheta/(3+lTheta), m2 = lPhi/(3+lTheta), m3 = lThetaPhi/(3+lTheta)
	    for(unsigned int ix=0;ix<3;++ix) {
	      val[ix] = moment[ix]->bins()[ibin].numEntries()>0 && moment[ix]->bins()[ibin].effNumEntries()>0 ? moment[ix]->bins()[ibin].mean(2)   : 0.;
	      err[ix] = moment[ix]->bins()[ibin].numEntries()>1 && moment[ix]->bins()[ibin].effNumEntries()>1 ? moment[ix]->bins()[ibin].stdErr(2) : 0.;
	    }
	    // values of the lambdas and their errors
	    double l1 = 3.*val[0]/(1.-val[0]);
	    double e1 = 3./sqr(1.-val[0])*err[0];
	    double l2 = (3.+l1)*val[1];
	    double e2 = 3./sqr(1.-val[0])*sqrt(sqr(err[0]*val[1])+sqr(err[1]*(1.-val[0])));
	    double l3 = (3.+l1)*val[2];
	    double e3 = 3./sqr(1.-val[0])*sqrt(sqr(err[0]*val[1])+sqr(err[1]*(1.-val[0])));
	    // fill the scatters
	    lTheta   ->bin(ibin).setVal(l1);
            lTheta   ->bin(ibin).setErr(make_pair(e1,e1));
	    lPhi     ->bin(ibin).setVal(l2);
	    lPhi     ->bin(ibin).setErr(make_pair(e2,e2));
	    lThetaPhi->bin(ibin).setVal(l3);
	    lThetaPhi->bin(ibin).setErr(make_pair(e3,e3));
	  }
       	}
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p_Jpsi[2][5][3];
    Histo1DGroupPtr _h_Jpsi;
    vector<double>  _ybins;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1244315);

}
