// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Charged-particle multiplicities, R, and event shapes between 12 and 43 GeV
  class TASSO_1984_I195333 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1984_I195333);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      declare(Sphericity(fs), "Sphericity");
      declare(Thrust(fs), "Thrust");

      // Counters for R
      book(_c_hadrons, "sigma_hadrons", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons,   "sigma_muons"  , refData<YODA::BinnedEstimate<string>>(1,1,1));
      for (const string& en : _c_hadrons.binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms[0] = en;
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[0] = en;
            break;
          }
        }
      }
      book(_h_weight, "/TMP/HWeight");
      unsigned int iloc(0);
      sqs = sqrtS();
      if (isCompatibleWithSqrtS(14*GeV)) {
        iloc = 1;
        sqs = 14.;
      }
      else if (isCompatibleWithSqrtS(22*GeV)) {
        iloc = 2;
        sqs = 22.;
      }
      else if (isCompatibleWithSqrtS(34*GeV)) {
        iloc = 3;
        sqs = 34.;
      }
      if (iloc!=0) {
        book(_h_mult, 3,1,iloc);
        book(_h_p,    5,1,iloc);
        book(_h_xp,   6,1,iloc);
        book(_h_pl,   7,1,iloc);
        book(_h_pt,   8,1,iloc);
        book(_h_pt2,  9,1,iloc);
        book(_h_xl,  10,1,iloc);
        book(_h_xT,  11,1,iloc);
        book(_h_S,   12,1,iloc);
        book(_h_T,   13,1,iloc);
        book(_h_y,   14,1,iloc);
      }
      // average event shapes
      book(_n_charged ,4 ,1,  1);
      book(_n_total   ,"TMP/ntotal",refData<YODA::BinnedEstimate<string>>(4 ,1,  2));
      book(_sphericity,4 ,1,  4);
      book(_thrust    ,4 ,1,  5);
      book(_p_total   ,4 ,1,  6);
      book(_p_l       ,4 ,1,  7);
      book(_pt        ,4 ,1,  8);
      book(_pt2       ,4 ,1,  9);
      book(_pt2_in    ,4 ,1, 10);
      book(_pt2_out   ,4 ,1, 11);
      for (const string& en : _n_charged.binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms[1] = en;
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[1] = en;
            break;
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      double ntotal(0);
      unsigned int nCharged(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        if ((p.pid()!=PID::PHOTON && p.abspid()!=PID::ELECTRON) ||
            p.parents().empty() || p.parents()[0].pid()!=PID::PI0) {
          ntotal += 1.;
        }
        else if (p.pid()==PID::PHOTON) {
          ntotal += 0.5 ;
        }
        else {
          ntotal += 0.25;
          nCharged -=1;
        }
        if (PID::isCharged(p.pid())) {
          ++nCharged;
        }
      }
      // mu+mu- + photons
      if (nCount[-13]==1 and nCount[13]==1 && int(fs.particles().size())==2+nCount[22]) {
        _c_muons->fill(_ecms[0]);
        return;
      }
      // everything else
      _c_hadrons->fill(_ecms[0]);
      _h_weight->fill();
      _n_charged->fill(_ecms[1],nCharged);
      _n_total->fill(_ecms[1],ntotal);
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _thrust->fill(_ecms[1],thrust.thrust());
      // sphericity
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _sphericity->fill(_ecms[1],sphericity.sphericity());
      // global distributions
      if (_h_mult)  _h_mult->fill(nCharged);
      if (_h_S)     _h_S   ->fill(sphericity.sphericity());
      if (_h_T)     _h_T   ->fill(thrust.thrust());
      // single particle distributions
      for (const Particle& p : fs.particles()) {
        if (!PID::isCharged(p.pid())) continue;
        const Vector3 mom3 = p.p3();
        double pp = mom3.mod();
        _p_total->fill(_ecms[1],pp);
        if (_h_p)  _h_p ->fill(pp);
        if (_h_xp) _h_xp->fill(2.*pp/sqrtS());
        const double mom = dot(sphericity.sphericityAxis(), mom3);
        _p_l->fill(_ecms[1],fabs(mom));
        if (_h_pl)  _h_pl->fill(fabs(mom));
        if (_h_xl)  _h_xl->fill(2.*fabs(mom)/sqrtS());
        const double pTin = dot(mom3, sphericity.sphericityMajorAxis());
        _pt2_in->fill(_ecms[1],sqr(pTin));
        const double pTout = dot(mom3, sphericity.sphericityMinorAxis());
        _pt2_out->fill(_ecms[1],sqr(pTout));
        double pT = sqr(pTin) + sqr(pTout);
        _pt2->fill(_ecms[1],pT);
        if (_h_pt2)  _h_pt2->fill(pT);
        pT=sqrt(pT);
        _pt->fill(_ecms[1],pT);
        if (_h_pt)  _h_pt->fill(pT);
        if (_h_xT)  _h_xT->fill(2.*pT/sqrtS());
        if (_h_y) {
          const double rap = 0.5 * log((p.E() + mom) / (p.E() - mom));
          _h_y->fill(fabs(rap));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 1);
      divide(_c_hadrons,_c_muons,mult);
      // charged particle multiplicity distribution
      if (_h_mult)  normalize(_h_mult,1.);
      // charged fraction
      book(mult,4,1,2);
      divide(_n_charged,_n_total,mult);
      book(mult,4,1,3);
      const auto & b = _n_charged->binAt(_ecms[1]);
      mult->binAt(_ecms[1]).set(b.stdDev(2),b.stdErr(2));
      // scale the distributions
      scale(_h_p  ,1./_h_weight->sumW());
      scale(_h_xp ,1./_h_weight->sumW());
      scale(_h_pl ,1./_h_weight->sumW());
      scale(_h_pt ,1./_h_weight->sumW());
      scale(_h_pt2,1./_h_weight->sumW());
      scale(_h_xl ,1./_h_weight->sumW());
      scale(_h_xT ,1./_h_weight->sumW());
      scale(_h_S  ,1./_h_weight->sumW());
      scale(_h_T  ,1./_h_weight->sumW());
      scale(_h_y  ,1./_h_weight->sumW());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_p,_h_xp,_h_pl,_h_pt,_h_pt2,_h_xl,_h_xT,_h_S,_h_T,_h_y;
    BinnedHistoPtr<int> _h_mult;
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    BinnedProfilePtr<string> _n_charged, _n_total,_sphericity, _thrust, _p_total, _p_l, _pt, _pt2, _pt2_in, _pt2_out;
    CounterPtr  _h_weight;
    double sqs;
    string _ecms[2];
    /// @}

  };



  RIVET_DECLARE_PLUGIN(TASSO_1984_I195333);

}
