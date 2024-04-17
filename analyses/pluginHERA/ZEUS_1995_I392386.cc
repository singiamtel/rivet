// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"


namespace Rivet {

const vector<double> QEdges {10., 20., 40., 80., 160., 320.};
const vector<double> xEdges {0.6e-3,1.2e-3,2.4e-3,1.0e-2,5.0e-2};

  /// @brief Charged particle multiplicity and momentum spectra in Breit frame at ZEUS
  class ZEUS_1995_I392386 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1995_I392386);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
        declare(DISKinematics(), "Kinematics");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
        const ChargedFinalState fs;
        declare(fs, "FS");


      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)

        for(int iQ = 0; iQ < 11; ++iQ) {
          book(_Nevt_after_cuts_Q[iQ], "TMP/Nevt_after_cuts_Q"+ to_string(iQ));
        }

        book(_p["mult1"], 1,1,1);
        book(_p["mult2"], 2,1,1);
        book(_p["mult3"], 3,1,1);
        book(_p["mult4"], 4,1,1);
        book(_p["mom1"], 5,1,1);
        book(_p["mom2"], 6,1,1);
        book(_p["mom3"], 7,1,1);
        book(_p["mom4"], 8,1,1);

        book(_i["nch1"], 9, 1, 1);  // Multiplicity
        book(_i["nch2"], 10, 1, 1);
        book(_i["nch3"], 10, 1, 2);
        book(_i["nch4"], 10, 1, 3);
        book(_i["nch5"], 11, 1, 1);
        book(_i["nch6"], 11, 1, 2);
        book(_i["nch7"], 11, 1, 3);
        book(_i["nch8"], 11, 1, 4);
        book(_i["nch9"], 12, 1, 1);
        book(_i["nch10"], 12, 1, 2);

        book(_h["loginvmom1"], 13, 1, 1); // Momentum spectra
        book(_h["loginvmom2"], 14, 1, 1);
        book(_h["loginvmom3"], 14, 1, 2);
        book(_h["loginvmom4"], 14, 1, 3);
        book(_h["loginvmom5"], 15, 1, 1);
        book(_h["loginvmom6"], 15, 1, 2);
        book(_h["loginvmom7"], 15, 1, 3);
        book(_h["loginvmom8"], 15, 1, 4);
        book(_h["loginvmom9"], 16, 1, 1);
        book(_h["loginvmom10"], 16, 1, 2);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "FS");
        const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");

        double xbj = dk.x(); // momentum fraction
        double Q2 = dk.Q2()/GeV; // momentum transfer
        const LorentzTransform Breitboost = dk.boostBreit();

        // Multiplicity counters
        int n911(0), n1011(0), n1012(0), n1013(0), n1111(0), n1112(0), n1113(0), n1114(0), n1211(0), n1212(0);


        if(0.6e-3<xbj && xbj<1.2e-3) {
            if(10<Q2 && Q2<20) {
                    _Nevt_after_cuts_Q[1] -> fill();
            }
            if(10<Q2 && Q2<20) {
                  _Nevt_after_cuts_Q[2] -> fill();
            }
            if(20<Q2 && Q2<40) {
                  _Nevt_after_cuts_Q[3] -> fill();
            }
            if(40<Q2 && Q2<80) {
                  _Nevt_after_cuts_Q[4] -> fill();
            }
        }
        if(2.4e-3<xbj && xbj<1.0e-2) {
            if(20<Q2 && Q2<40) {
                  _Nevt_after_cuts_Q[5] -> fill();
            }
            if(40<Q2 && Q2<80){
                  _Nevt_after_cuts_Q[6] -> fill();
            }
            if(80<Q2 && Q2<160) {
                  _Nevt_after_cuts_Q[7] -> fill();
            }
            if(160<Q2 && Q2<320) {
                 _Nevt_after_cuts_Q[8] -> fill();
            }
        }
        if(1.0e-2<xbj && xbj<5.0e-2) {
            if(320<Q2 && Q2<640) {
                 _Nevt_after_cuts_Q[9] -> fill();
            }
            if(640<Q2 && Q2<1280) {
                 _Nevt_after_cuts_Q[10] -> fill();
            }
        }


        for (const Particle& p : cfs.particles()) {
            //??? calculating ln(1/x_p) ??? --> part to ask
            const FourMomentum BrMom = Breitboost.transform(p.momentum());
            double pp = sqrt(BrMom.px()*BrMom.px() + BrMom.py()*BrMom.py() + BrMom.pz()*BrMom.pz() );
            double xp = 2*pp/(sqrt(Q2));
            const double logInvScaledMom = log(1/xp);

            if ( BrMom.pz() > 0. ) continue;

            if(0.6e-3<xbj && xbj<1.2e-3) {
                _p["mom1"] ->fill(Q2, logInvScaledMom);
                if(10<Q2 && Q2<20) {
                    _h["loginvmom1"] ->fill(logInvScaledMom);
                    ++n911;
                }
            }

          if(1.2e-3<xbj && xbj<2.4e-3) {
                _p["mom2"] ->fill(Q2, logInvScaledMom);
              if(10<Q2 && Q2<20) {
                  _h["loginvmom2"] ->fill(logInvScaledMom);
                  ++n1011;
                  }
              if(20<Q2 && Q2<40) {
                  _h["loginvmom3"] ->fill(logInvScaledMom);
                  ++n1012;
                  }
              if(40<Q2 && Q2<80) {
                  _h["loginvmom4"] ->fill(logInvScaledMom);
                  ++n1013;
                  }
            }

            if(2.4e-3<xbj && xbj<1.0e-2) {
                _p["mom3"] ->fill(Q2, logInvScaledMom);
                if(20<Q2 && Q2<40) {
                    _h["loginvmom5"] ->fill(logInvScaledMom);
                    ++n1111;
                    }
                if(40<Q2 && Q2<80){
                    _h["loginvmom6"] ->fill(logInvScaledMom);
                    ++n1112;
                    }
                if(80<Q2 && Q2<160) {
                    _h["loginvmom7"] ->fill(logInvScaledMom);
                    ++n1113;
                    }
                if(160<Q2 && Q2<320) {
                    _h["loginvmom8"] ->fill(logInvScaledMom);
                    ++n1114;
                    }
            }

            if(1.0e-2<xbj && xbj<5.0e-2) {
                _p["mom4"] ->fill(Q2,logInvScaledMom);
                if(320<Q2 && Q2<640) {
                    _h["loginvmom9"] ->fill(logInvScaledMom);
                    ++n1211;
                    }
                    if(640<Q2 && Q2<1280) {
                    _h["loginvmom10"] ->fill(logInvScaledMom);
                    ++n1212;
                    }
            }

            }

        if(0.6e-3<xbj && xbj<1.2e-3) {
            if(10<Q2 && Q2<20) {
                _p["mult1"] ->fill(Q2, n911);
                _i["nch1"] ->fill(n911);
            }
        }

      if(1.2e-3<xbj && xbj<2.4e-3) {
          if(10<Q2 && Q2<80) {
            _p["mult2"] ->fill(Q2, n1011+n1012+n1013);
              if(10<Q2 && Q2<20) {
                  _i["nch2"] ->fill(n1011); }
              if(20<Q2 && Q2<40) {
                  _i["nch3"] ->fill(n1012);}
              if(40<Q2 && Q2<80) {
                  _i["nch4"] ->fill(n1013);}
          }
        }

        if(2.4e-3<xbj && xbj<1.0e-2) {
            _p["mult3"] ->fill(Q2,n1111+n1112+n1113+n1114);
            if(20<Q2 && Q2<40) {
                _i["nch5"] ->fill(n1111);}
            if(40<Q2 && Q2<80) {
                _i["nch6"] ->fill(n1112);}
            if(80<Q2 && Q2<160) {
                _i["nch7"] ->fill(n1113);}
            if(160<Q2 && Q2<320) {
                _i["nch8"] ->fill(n1114);}

        }

        if(1.0e-2<xbj && xbj<5.0e-2) {
            _p["mult4"] ->fill(Q2, n1211+n1212);
            if(320<Q2 && Q2<640) {
                _i["nch9"] ->fill(n1211);}
            if(640<Q2 && Q2<1280) {
                _i["nch10"] ->fill(n1212); }
        }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
        normalize(_i["nch1"]); //multiplicity
        normalize(_i["nch2"]);
        normalize(_i["nch3"]);
        normalize(_i["nch4"]);
        normalize(_i["nch5"]);
        normalize(_i["nch6"]);
        normalize(_i["nch7"]);
        normalize(_i["nch8"]);
        normalize(_i["nch9"]);
        normalize(_i["nch10"]);

        if(dbl(*_Nevt_after_cuts_Q[1])>0 ) scale(_h["loginvmom1"],1./ *_Nevt_after_cuts_Q[1]);  //momentum
        if(dbl(*_Nevt_after_cuts_Q[2])>0 ) scale(_h["loginvmom2"],1./ *_Nevt_after_cuts_Q[2]);
        if(dbl(*_Nevt_after_cuts_Q[3])>0 ) scale(_h["loginvmom3"],1./ *_Nevt_after_cuts_Q[3] );
        if(dbl(*_Nevt_after_cuts_Q[4])>0 ) scale(_h["loginvmom4"],1./ *_Nevt_after_cuts_Q[4] );
        if(dbl(*_Nevt_after_cuts_Q[5])>0 ) scale(_h["loginvmom5"],1./ *_Nevt_after_cuts_Q[5] );
        if(dbl(*_Nevt_after_cuts_Q[6])>0 ) scale(_h["loginvmom6"],1./ *_Nevt_after_cuts_Q[6]);
        if(dbl(*_Nevt_after_cuts_Q[7])>0 ) scale(_h["loginvmom7"],1./ *_Nevt_after_cuts_Q[7]);
        if(dbl(*_Nevt_after_cuts_Q[8])>0 ) scale(_h["loginvmom8"],1./ *_Nevt_after_cuts_Q[8]);
        if(dbl(*_Nevt_after_cuts_Q[9])>0 ) scale(_h["loginvmom9"],1./ *_Nevt_after_cuts_Q[9] );
        if(dbl(*_Nevt_after_cuts_Q[10])>0 ) scale(_h["loginvmom10"],1./ *_Nevt_after_cuts_Q[10] );

    }

    ///@}


    /// @name Histograms
    ///@{

    Estimate1DPtr _h_mult1;
    Estimate1DPtr _h_mult2;
    Estimate1DPtr _h_mult3;
    Estimate1DPtr _h_mult4;
    Estimate1DPtr _h_mom1;
    Estimate1DPtr _h_mom2;
    Estimate1DPtr _h_mom3;
    Estimate1DPtr _h_mom4;

    map<string, Histo1DPtr> _h;
    map<string, BinnedHistoPtr<int> > _i;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    CounterPtr _Nevt_after_cuts_Q[11];

    ///@}


  };


  RIVET_DECLARE_PLUGIN(ZEUS_1995_I392386);

}
