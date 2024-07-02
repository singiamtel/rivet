// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"


namespace Rivet {


  /// @brief  Lambda_c+ -> Lambda0 l+ nu_l
  class BESIII_2023_I2689064 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2689064);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::K0S);
      LAMBDAC.addStable(PID::ETA);
      LAMBDAC.addStable(PID::ETAPRIME);
      declare(LAMBDAC, "LAMBDAC");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_q[ix], 1, 1+ix, 1);
        for (unsigned int iy=0; iy<4; ++iy) {
          book(_h_kin[ix][iy], 3, ix+1, iy+1);
          if (iy<2) {
            book(_h_asym_l[ix][iy],"/TMP/h_asym_l_"+toString(ix+1)+"_"+toString(iy+1),refData(1,ix+1,2));
            book(_h_asym_p[ix][iy],"/TMP/h_asym_p_"+toString(ix+1)+"_"+toString(iy+1),refData(2,1,ix+1));
          }
        }
      }
      book(_b_ctheta, {0.,0.2,0.4,0.6,0.8,1.0,1.2,1.37});
      for (auto& b : _b_ctheta->bins()) {
        book(b, "/TMP/hLambda_"+toString(b.index()), 20, -1., 1.);
      }
      book(_c,"/TMP/n_lambda");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for (unsigned int ix=0; ix<LAMBDAC.decaying().size(); ++ix) {
        _c->fill();
        unsigned int il=0;
        if      (LAMBDAC.modeMatches(ix,4,mode1) ) il=0;
        else if (LAMBDAC.modeMatches(ix,4,mode2) ) il=1;
        else {
          continue;
        }
      	const Particle& pp = LAMBDAC.decayProducts()[ix].at(2212)[0];
      	const Particle& pim= LAMBDAC.decayProducts()[ix].at(-211)[0];
      	const Particle& ep = LAMBDAC.decayProducts()[ix].at( -11-2*il)[0];
      	const Particle& nue= LAMBDAC.decayProducts()[ix].at(  12+2*il)[0];
      	if (LAMBDAC.decaying()[ix].children(Cuts::pid==PID::LAMBDA).empty()) continue;
      	FourMomentum pLambda = pp.mom()+pim.mom();
       	FourMomentum qq = LAMBDAC.decaying()[ix].mom()-pLambda;
        double q2 = qq.mass2();
        _h_kin[il][0]->fill(q2);
        _h_q[il]->fill(q2);
        // boost momenta to LAMBDAC rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(LAMBDAC.decaying()[ix].mom().betaVec());
        FourMomentum pLam = boost.transform(pLambda);
        Matrix3 ptoz(-pLam.p3().unit(), Vector3(0,0,1));
        boost.preMult(ptoz);
        // the momenta in frane to W along z
        FourMomentum pD  = boost.transform(LAMBDAC.decaying()[ix].mom());
        FourMomentum pP  = boost.transform(pp .mom());
        FourMomentum ppi = boost.transform(pim.mom());
        FourMomentum pe  = boost.transform(ep .mom());
        FourMomentum pnu = boost.transform(nue.mom());
        pLambda = pP+ppi;
        qq = pD-pLambda;
        LorentzTransform boostL = LorentzTransform::mkFrameTransformFromBeta(pLambda.betaVec());
        Vector3 axisP = boostL.transform(pP).p3().unit();
        const double cThetaP = axisP.dot(pLambda.p3().unit());
        _h_kin[il][1]->fill(cThetaP);
        _b_ctheta->fill(q2,cThetaP);
        if (cThetaP>0.) _h_asym_p[il][0]->fill(q2);
        else            _h_asym_p[il][1]->fill(q2);
        LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
        Vector3 axisE = boostW.transform(pe).p3().unit();
        double cThetaE = -axisE.dot(qq.p3().unit());
        _h_kin[il][2]->fill(cThetaE);
        if (cThetaE>0.) _h_asym_l[il][0]->fill(q2);
        else            _h_asym_l[il][1]->fill(q2);
        axisP.setZ(0.);
        axisE.setZ(0.);
        double chi = atan2(axisE.cross(axisP).dot(qq.p3().unit()), axisE.dot(axisP));
        _h_kin[il][3]->fill(chi);
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if (hist->numEntries()==0.)  return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (const auto& bin : hist->bins()) {
        const double Oi = bin.sumW();
        if (Oi==0.) continue;
        const double ai = 0.5*(bin.xMax()-bin.xMin());
        const double bi = 0.5*ai*(bin.xMax()+bin.xMin());
        const double Ei = bin.errW();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double tau = 201.5e-3; // lifetime in ps from PDG 2023
      const double br = 0.641; // br for lambda -> p pi from PDG 2023
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_q[ix], 1.0/tau/br/ *_c);
        normalize(_h_kin[ix], 1.0, false);
      }
      // asymmetries
      Estimate1DPtr ap[2];
      for (unsigned int ix=0; ix<2; ++ix) {
        Estimate1DPtr as;
        book(as, 1, 1+ix, 2);
        asymm(_h_asym_l[ix][1], _h_asym_l[ix][0], as);
        book(ap[ix], 2, 1, 1+ix);
        asymm(_h_asym_p[ix][0], _h_asym_p[ix][1], ap[ix]);
      }
      Estimate1DPtr a2;
      book(a2, 1, 3, 2);
      divide(ap[1], ap[0], a2);
      // ratios
      Estimate1DPtr r1;
      book(r1,1,3,1);
      divide(_h_q[1], _h_q[0],r1);
      // alpha parameters
      Estimate1DPtr _h_alpha;
      book(_h_alpha,2,2,1);
      double lambda = 0.748; // PDG 2023
      normalize(_b_ctheta);
      for(unsigned int ix=0;ix<7;++ix) {
	pair<double,double> alpha = calcAlpha(_b_ctheta->bin(ix+1));
	alpha.first  /= lambda;
	alpha.second /= lambda;
	_h_alpha->bin(ix+1).set(alpha.first, alpha.second);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_kin[2][4], _h_q[2] ,_h_asym_l[2][2], _h_asym_p[2][2];
    Histo1DGroupPtr _b_ctheta;
    CounterPtr _c;
    const map<PdgId,unsigned int> mode1 = { { 2212,1}, {-211,1}, {-11,1}, { 12,1}};
    const map<PdgId,unsigned int> mode2 = { { 2212,1}, {-211,1}, {-13,1}, { 14,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2689064);

}
