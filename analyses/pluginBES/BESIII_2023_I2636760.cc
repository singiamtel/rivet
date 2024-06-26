// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Sigma + -> p gamma
  class BESIII_2023_I2636760 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2636760);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");

      // Book histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_n[ix],"TMP/n_" + toString(ix+1));
        book(_t[ix],"TMP/t_" + toString(ix+1));
      }
      book(_n[2],"TMP/n_3");
      book(_t[2],"TMP/t_3");
    }


    void findChildren(const Particle& p,map<long,int>& nRes, int& ncount) {
      for (const Particle& child : p.children()) {
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
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.pid()>0) {
        axis = beams.first .mom().p3().unit();
      }
      else {
	      axis = beams.second.mom().p3().unit();
      }
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // loop over Sigma+ baryons
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      Particle Sigma, SigBar;
      bool matched(false);
      for (const Particle& p :  ufs.particles(Cuts::abspid==3222)) {
       	if (p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p, nRes, ncount);
       	matched=false;
       	// check for antiparticle
      	for (const Particle& p2 :  ufs.particles(Cuts::pid==-p.pid())) {
      	  if(p2.children().empty()) continue;
      	  map<long,int> nRes2=nRes;
      	  int ncount2 = ncount;
      	  findChildren(p2, nRes2, ncount2);
      	  if (ncount2==0) {
      	    matched = true;
      	    for (const auto& val : nRes2) {
      	      if (val.second!=0) {
                matched = false;
                break;
      	      }
      	    }
      	    // fond baryon and antibaryon
      	    if (matched) {
              if (p.pid()>0) {
                Sigma = p;
                SigBar = p2;
              }
              else {
                Sigma = p2;
                SigBar = p;
              }
       	      break;
       	    }
       	  }
       	}
      	if (matched) break;
      }
      if (!matched) vetoEvent;
      if (Sigma .children().size()!=2) vetoEvent;
      if (SigBar.children().size()!=2) vetoEvent;
      // find proton
      bool radiative[2]={false,false};
      // identify Sigma decay
      Particle baryon1;
      if (Sigma.children()[0].pid()==PID::PROTON && Sigma.children()[1].pid()==PID::PI0) {
        radiative[0]=false;
        baryon1 = Sigma.children()[0];
      }
      else if (Sigma.children()[1].pid()==PID::PROTON && Sigma.children()[0].pid()==PID::PI0) {
        radiative[0]=false;
        baryon1 = Sigma.children()[1];
      }
      else if (Sigma.children()[0].pid()==PID::PROTON && Sigma.children()[1].pid()==PID::PHOTON) {
        radiative[0]=true;
        baryon1 = Sigma.children()[0];
      }
      else if (Sigma.children()[1].pid()==PID::PROTON && Sigma.children()[0].pid()==PID::PHOTON ) {
        radiative[0]=true;
        baryon1 = Sigma.children()[1];
      }
      else {
        vetoEvent;
      }
      // antisigma decay
      Particle baryon2;
      if (SigBar.children()[0].pid()==PID::ANTIPROTON && SigBar.children()[1].pid()==PID::PI0 ) {
        radiative[1]=false;
        baryon2 = SigBar.children()[0];
      }
      else if (SigBar.children()[1].pid()==PID::ANTIPROTON && SigBar.children()[0].pid()==PID::PI0 ) {
        radiative[1]=false;
        baryon2 = SigBar.children()[1];
      }
      else if (SigBar.children()[0].pid()==PID::ANTIPROTON && SigBar.children()[1].pid()==PID::PHOTON ) {
        radiative[1]=true;
        baryon2 = SigBar.children()[0];
      }
      else if (SigBar.children()[1].pid()==PID::ANTIPROTON && SigBar.children()[0].pid()==PID::PHOTON ) {
        radiative[1]=true;
        baryon2 = SigBar.children()[1];
      }
      else {
        vetoEvent;
      }
      if (radiative[0] == radiative[1]) vetoEvent;
      // boost to the Sigma rest frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Sigma.mom().betaVec());
      Vector3 e1z = Sigma.mom().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      Vector3 axis1 = boost1.transform(baryon1.mom()).p3().unit();
      double n1x(e1x.dot(axis1)), n1z(e1z.dot(axis1));
      // boost to the Sigma bar
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(SigBar.mom().betaVec());
      Vector3 axis2 = boost2.transform(baryon2.mom()).p3().unit();
      double n2x(e1x.dot(axis2)),n2z(e1z.dot(axis2));
      double cosL = axis.dot(Sigma.mom().p3().unit());
      double sinL = sqrt(1.-sqr(cosL));
      double T1 = sqr(sinL)*n1x*n2x+sqr(cosL)*n1z*n2z;
      // sigma+ -> p gamma
      if (radiative[0]) {
        _n[0]->fill();
        _n[2]->fill();
        _t[0]->fill(T1);
        _t[2]->fill(T1);
      }
      // sigmabar- -> pbar gamma
      else {
        _n[1]->fill();
        _n[2]->fill();
        _t[1]->fill(T1);
        _t[2]->fill(T1);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // values of constants
      const double aPsi  =-0.508;
      const double aPlus =-0.998;
      const double factor = 45.*(3. +aPsi)/(11. + 5.*aPsi)/aPlus;
      // alpha from the moments
      for (unsigned int ix=0;ix<3;++ix) {
        double value = _t[ix]->val()/_n[ix]->val();
        double error = _t[ix]->err()/_n[ix]->val();
        value *= factor;
        error *= abs(factor);
        if (ix==1) value *=-1.;
        Estimate0DPtr alpha;
        book(alpha, 1, 1, 1+ix);
        alpha->set(value, error);
      }

    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _n[3];
    Histo1DPtr _h_ctheta[3];
    CounterPtr _t[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2636760);

}
