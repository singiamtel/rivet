// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> D* semileptonic
  class BELLE_2023_I2624324 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2023_I2624324);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521), "UFS");
      // histograms
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h_aver[ix], 2, 1+ix, 1);
        for (unsigned int iy=0; iy<4; ++iy) {
          book(_h[ix][iy], 1, 1+ix, 1+iy);
        }
      }
    }

    /// Perform the per-event analysis
    bool analyzeDecay(const Particle& mother, const vector<int>& ids) {
      // There is no point in looking for decays with less particles than to be analysed
      if (mother.children().size() == ids.size()) {
        bool decayfound = true;
        for (int id : ids) {
          if (!contains(mother, id)) decayfound = false;
        }
        return decayfound;
      }
      return false;
    }

    bool contains(const Particle& mother, int id) {
      return any(mother.children(), HasPID(id));
    }

    double recoilW(const Particle& mother) {
      FourMomentum lepton, neutrino, meson, q;
      for(const Particle& c : mother.children()) {
        if (c.isNeutrino()) neutrino=c.mom();
        if (c.isLepton() && !c.isNeutrino()) lepton =c.mom();
        if (c.isHadron()) meson=c.mom();
      }
      q = lepton + neutrino; //no hadron before
      double mb2= mother.mom()*mother.mom();
      double mD2 = meson*meson;
      return (mb2 + mD2 - q*q )/ (2. * sqrt(mb2) * sqrt(mD2) );
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      FourMomentum pl, pnu, pB, pD, pDs, ppi;
      // Iterate of B mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        pB = p.mom();
        // Find semileptonic decays
        int sign = p.pid()/p.abspid();
        int iDStar = sign*(p.abspid()==511 ? -413 : -423);
        int iloc=-1;
        if (analyzeDecay(p, {iDStar,12*sign,-11*sign})) iloc = 0;
        else if(analyzeDecay(p, {iDStar,14*sign,-13*sign})) iloc = 1;
        else continue;
        if (p.abspid()==521) iloc+=2;
        double w = recoilW(p);
        _h[0][iloc]->fill(w);
        _h_aver[0] ->fill(w);
        // Get the necessary momenta for the angles
        bool foundDdecay=false;
        for (const Particle& c : p.children()) {
          if (c.abspid()==413 || c.abspid()==423) {
            if ((c.pid() == -413  && (analyzeDecay(c, {PID::PIMINUS, PID::D0BAR}) || analyzeDecay(c, {PID::PI0, PID::DMINUS})) ) ||
                (c.pid() ==  413  && (analyzeDecay(c, {PID::PIPLUS , PID::D0   }) || analyzeDecay(c, {PID::PI0, PID::DPLUS })) ) ||
                (c.pid() == -423  && analyzeDecay(c, {PID::PI0, PID::D0BAR }) ) ||
                (c.pid() ==  423  && analyzeDecay(c, {PID::PI0, PID::D0    }) )) {
              foundDdecay=true;
              pDs = c.mom();
              for (const Particle & dc : c.children()) {
                if (dc.hasCharm()) pD = dc.mom();
                else ppi = dc.mom();
              }
            }
          }
          else if (c.abspid() ==  11 || c.abspid() ==  13) pl  = c.mom();
          else if (c.abspid() ==  12 || c.abspid() ==  14) pnu = c.mom();
        }
        // This is the angle analysis
        if (!foundDdecay) continue;
        // First boost all relevant momenta into the B-rest frame
        const LorentzTransform B_boost = LorentzTransform::mkFrameTransformFromBeta(pB.betaVec());
        // Momenta in B rest frame:
        FourMomentum lv_brest_Dstar = B_boost.transform(pDs);
        FourMomentum lv_brest_w     = B_boost.transform(pB - pDs);
        FourMomentum lv_brest_D     = B_boost.transform(pD);
        FourMomentum lv_brest_lep   = B_boost.transform(pl);

        const LorentzTransform Ds_boost = LorentzTransform::mkFrameTransformFromBeta(lv_brest_Dstar.betaVec());
        FourMomentum lv_Dstarrest_D     = Ds_boost.transform(lv_brest_D);
        const LorentzTransform W_boost  = LorentzTransform::mkFrameTransformFromBeta(lv_brest_w.betaVec());
        FourMomentum lv_wrest_lep       = W_boost.transform(lv_brest_lep);

        double cos_thetaV = cos(lv_brest_Dstar.p3().angle(lv_Dstarrest_D.p3()));
        _h[2][iloc]->fill(cos_thetaV);
        _h_aver[2] ->fill(cos_thetaV);

        double cos_thetaL = cos(lv_brest_w.p3().angle(lv_wrest_lep.p3()));
        _h[1][iloc]->fill(cos_thetaL);
        _h_aver[1] ->fill(cos_thetaL);

        Vector3 LTrans = lv_wrest_lep.p3()   - cos_thetaL*lv_wrest_lep.p3().perp()*lv_brest_w.p3().unit();
        Vector3 VTrans = lv_Dstarrest_D.p3() - cos_thetaV*lv_Dstarrest_D.p3().perp()*lv_brest_Dstar.p3().unit();
        float chi = atan2(LTrans.cross(VTrans).dot(lv_brest_w.p3().unit()), LTrans.dot(VTrans));
        if (chi<0.) chi+=2.*M_PI;
      	_h[3][iloc]->fill(chi);
      	_h_aver[3] ->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h     , 1.0);
      normalize(_h_aver, 1.0);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][4];
    Histo1DPtr _h_aver[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2023_I2624324);

}
