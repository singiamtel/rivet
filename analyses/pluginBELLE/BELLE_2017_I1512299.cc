// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> D*+ semileptonic
  class BELLE_2017_I1512299 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1512299);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_w     , "TMP/h_w"     , refData(1, 1, 1));
      book(_h_costhv, "TMP/h_costhv", refData(2, 1, 1));
      book(_h_costhl, "TMP/h_costhl", refData(3, 1, 1));
      book(_h_chi   , "TMP/h_chi"   , refData(4, 1, 1));
      book(_nB,"TMP/nB");
    }


    /// Perform the per-event analysis
    bool analyzeDecay(Particle mother, vector<int> ids) {
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

    bool contains(Particle& mother, int id) {
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
      // Iterate of B0bar mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==-511)) {
        if(p.children().size()==1 && p.children()[0].abspid()==511) continue;
        _nB->fill();
        pB = p.momentum();
        // Find semileptonic decays
        if (analyzeDecay(p, {PID::DSTARPLUS,-12,11}) || analyzeDecay(p, {PID::DSTARPLUS,-14,13}) ) {
          _h_w->fill(recoilW(p));
          // Get the necessary momenta for the angles
          bool foundDdecay=false;
          for (const Particle & c : p.children()) {
            if ( (c.pid() == PID::DSTARPLUS)  && (analyzeDecay(c, {PID::PIPLUS, PID::D0}) || analyzeDecay(c, {PID::PI0, PID::DPLUS})) ) {
              foundDdecay=true;
              pDs = c.momentum();
              for (const Particle & dc : c.children()) {
                if (dc.hasCharm()) pD = dc.momentum();
                else ppi = dc.momentum();
              }
            }
            if (c.pid() ==  11 || c.pid() ==  13) pl  = c.momentum();
            if (c.pid() == -12 || c.pid() == -14) pnu = c.momentum();
          }
          // This is the angle analysis
          if (foundDdecay) {

            // First boost all relevant momenta into the B-rest frame
            const LorentzTransform B_boost = LorentzTransform::mkFrameTransformFromBeta(pB.betaVec());
            // Momenta in B rest frame
            FourMomentum lv_brest_Dstar = B_boost.transform(pDs);
            FourMomentum lv_brest_w     = B_boost.transform(pB - pDs);
            FourMomentum lv_brest_D     = B_boost.transform(pD);
            FourMomentum lv_brest_lep   = B_boost.transform(pl);
            
            const LorentzTransform Ds_boost = LorentzTransform::mkFrameTransformFromBeta(lv_brest_Dstar.betaVec());
            FourMomentum lv_Dstarrest_D     = Ds_boost.transform(lv_brest_D);
            const LorentzTransform W_boost  = LorentzTransform::mkFrameTransformFromBeta(lv_brest_w.betaVec());
            FourMomentum lv_wrest_lep       = W_boost.transform(lv_brest_lep);

            double cos_thetaV = cos(lv_brest_Dstar.p3().angle(lv_Dstarrest_D.p3()));
            _h_costhv->fill(cos_thetaV);

            double cos_thetaL = cos(lv_brest_w.p3().angle(lv_wrest_lep.p3()));
            _h_costhl->fill(cos_thetaL);

            Vector3 LTrans = lv_wrest_lep.p3()   - cos_thetaL*lv_wrest_lep.p3().perp()*lv_brest_w.p3().unit();
            Vector3 VTrans = lv_Dstarrest_D.p3() - cos_thetaV*lv_Dstarrest_D.p3().perp()*lv_brest_Dstar.p3().unit();
            float chi = atan2(LTrans.cross(VTrans).dot(lv_brest_w.p3().unit()), LTrans.dot(VTrans));
            if(chi < 0) chi += TWOPI;

            _h_chi->fill(chi);
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double GAMMA_B0 = 4.32e-13 *1e15; // Normalise histos to partial width ifetime converted to GeV (x10^-15)
      // factor 1/2 as e+mu in hists 
      scale(_h_w,      0.5*GAMMA_B0/ *_nB);
      scale(_h_costhv, 0.5*GAMMA_B0/ *_nB);
      scale(_h_costhl, 0.5*GAMMA_B0/ *_nB);
      scale(_h_chi,    0.5*GAMMA_B0/ *_nB);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_h_w,tmp);
      book(tmp,2,1,1);
      barchart(_h_costhv,tmp);
      book(tmp,3,1,1);
      barchart(_h_costhl,tmp);
      book(tmp,4,1,1);
      barchart(_h_chi,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_w;
    Histo1DPtr _h_costhv;
    Histo1DPtr _h_costhl;
    Histo1DPtr _h_chi;
    CounterPtr _nB;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2017_I1512299);


}
