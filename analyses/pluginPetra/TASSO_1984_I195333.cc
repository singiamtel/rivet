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
      book(_c_hadrons, "sigma_hadrons");
      book(_c_muons, "sigma_muons");
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
        _c_muons->fill();
        return;
      }
      // everything else
      _c_hadrons->fill();
      _h_weight->fill();
      _n_charged.fill(nCharged);
      _n_total.fill(ntotal);
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _thrust.fill(thrust.thrust());
      // sphericity
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _sphericity.fill(sphericity.sphericity());
      // global distributions
      if (_h_mult)  _h_mult->fill(nCharged);
      if (_h_S)     _h_S   ->fill(sphericity.sphericity());
      if (_h_T)     _h_T   ->fill(thrust.thrust());
      // single particle distributions
      for (const Particle& p : fs.particles()) {
        if (!PID::isCharged(p.pid())) continue;
        const Vector3 mom3 = p.p3();
        double pp = mom3.mod();
        _p_total.fill(pp);
        if (_h_p)  _h_p ->fill(pp);
        if (_h_xp) _h_xp->fill(2.*pp/sqrtS());
        const double mom = dot(sphericity.sphericityAxis(), mom3);
        _p_l.fill(fabs(mom));
        if (_h_pl)  _h_pl->fill(fabs(mom));
        if (_h_xl)  _h_xl->fill(2.*fabs(mom)/sqrtS());
        const double pTin = dot(mom3, sphericity.sphericityMajorAxis());
        _pt2_in.fill(sqr(pTin));
        const double pTout = dot(mom3, sphericity.sphericityMinorAxis());
        _pt2_out.fill(sqr(pTout));
        double pT = sqr(pTin) + sqr(pTout);
        _pt2.fill(pT);
        if (_h_pt2)  _h_pt2->fill(pT);
        pT=sqrt(pT);
        _pt.fill(pT);
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
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      scale({_c_hadrons, _c_muons}, fact);

      Estimate0D R = *_c_hadrons/ *_c_muons;
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      size_t idx = 0;
      if (isCompatibleWithSqrtS(12*GeV))          idx = 1;
      else if (isCompatibleWithSqrtS(14*GeV))     idx = 2;
      else if (isCompatibleWithSqrtS(22*GeV))     idx = 3;
      else if (isCompatibleWithSqrtS(25*GeV))     idx = 4;
      else if (inRange(sqrtS()/GeV, 27.4, 27.7))  idx = 5;
      else if (isCompatibleWithSqrtS(30.1*GeV))   idx = 6;
      else if (inRange(sqrtS()/GeV, 30.5, 31.5))  idx = 7;
      else if (inRange(sqrtS()/GeV, 32.5, 33.5))  idx = 8;
      else if (inRange(sqrtS()/GeV, 33.5, 34.5))  idx = 9;
      else if (inRange(sqrtS()/GeV, 34.5, 35.5))  idx = 10;
      else if (inRange(sqrtS()/GeV, 35.5, 36.7))  idx = 11;
      else if (inRange(sqrtS()/GeV, 38.7, 43.1))  idx = 12;
      mult->bin(idx).set(R.val(), R.errPos());

      // charged particle multiplicity distribution
      if (_h_mult)  normalize(_h_mult,2.);
      for (unsigned int iy=1;iy<12;++iy) {
        double value = 0.0, error = 0.0;
        if (iy==1) {
          value = _n_charged.xMean();
          error = _n_charged.xStdErr();
        }
        else if (iy==2) {
          double num = _n_charged.xMean();
          double den =   _n_total.xMean();
          value = num/den;
          error = value*sqrt(sqr(_n_charged.xStdErr()/num)+sqr(_n_total.xStdErr()/den));
        }
        else if (iy==3) {
          value = _n_charged.xStdDev();
          error = _n_charged.xStdErr();
        }
        else if (iy==4) {
          value = _sphericity.xMean();
          error = _sphericity.xStdErr();
        }
        else if (iy==5) {
          value = _thrust.xMean();
          error = _thrust.xStdErr();
        }
        else if (iy==6) {
          value = _p_total.xMean();
          error = _p_total.xStdErr();
        }
        else if (iy==7) {
          value = _p_l.xMean();
          error = _p_l.xStdErr();
        }
        else if (iy==8) {
          value = _pt.xMean();
          error = _pt.xStdErr();
        }
        else if (iy==9) {
          value = _pt2.xMean();
          error = _pt2.xStdErr();
        }
        else if (iy==10) {
          value = _pt2_in.xMean();
          error = _pt2_in.xStdErr();
        }
        else if (iy==11) {
          value = _pt2_out.xMean();
          error = _pt2_out.xStdErr();
        }
        BinnedEstimatePtr<string> mult;
        book(mult, 4, 1, iy);
        for (auto& b : mult->bins()) {
          const double Ecm = std::stod(b.xEdge());
          if (isCompatibleWithSqrtS(Ecm*GeV)) {
            b.set(value, error);
          }
        }
      }
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
    CounterPtr _c_hadrons, _c_muons;
    YODA::Dbn1D _n_charged, _n_total, _sphericity, _thrust, _p_total, _p_l, _pt, _pt2, _pt2_in, _pt2_out;
    CounterPtr  _h_weight;
    double sqs;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(TASSO_1984_I195333);

}
