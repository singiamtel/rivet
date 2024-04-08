// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief K+- Lambda asymmetries
  class DELPHI_1995_I382285 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1995_I382285);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_Kp, "/TMP/cos_Kp",20,-1.,1.);
      book(_h_Km, "/TMP/cos_Km",20,-1.,1.);
      book(_h_lm, "/TMP/cos_lm",20,-1.,1.);
      book(_h_lb, "/TMP/cos_lb",20,-1.,1.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      Vector3 beamAxis;
      if (beams.first.pid()==11) {
        beamAxis = beams.first .momentum().p3().unit();
      }
      else {
        beamAxis = beams.second.momentum().p3().unit();
      }
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle & p : ufs.particles(Cuts::abspid==3122 || Cuts::abspid==321 )) {
        double modp = p.momentum().p3().mod();
        if (p.abspid()==321) {
          if (modp<10. || modp>18.) continue;
          const double cK = beamAxis.dot(p.momentum().p3().unit());
          if (p.pid()>0)  _h_Kp->fill(cK);
          else            _h_Km->fill(cK);
        }
        else {
          if (modp<11.41 || modp>22.82) continue;
          const double cLam = beamAxis.dot(p.momentum().p3().unit());
          if (p.pid()>0)  _h_lm->fill(cLam);
          else            _h_lb->fill(cLam);
        }
      }
    }

    pair<double,double> calcAsymmetry(const Estimate1DPtr& hist) {
      double sum1(0.), sum2(0.);
      for (const auto& bin : hist->bins()) {
        double Oi = bin.val();
        if (Oi==0.) continue;
        const double bi = 4.*(bin.xMax()+bin.xMin())/(3.+sqr(bin.xMax())+bin.xMax()*bin.xMin()+sqr(bin.xMin()));
        const double Ei = bin.errAvg();
        sum1 += sqr(bi/Ei);
        sum2 += bi/sqr(Ei)*Oi;
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
       	normalize(_h_Kp);
       	normalize(_h_Km);
       	Estimate1DPtr sK;
        book(sK,"a_K",20,-1.,1.);
       	asymm(_h_Kp,_h_Km,sK);
        pair<double,double> alpha = calcAsymmetry(sK);
        BinnedEstimatePtr<string> h_K;
        book(h_K, 1,1,1);
       	h_K->bin(1).set(-alpha.first, alpha.second);

       	normalize(_h_lm);
       	normalize(_h_lb);
       	Estimate1DPtr sLam;
        book(sLam,"a_Lam",20,-1.,1.);
       	asymm(_h_lm,_h_lb,sLam);
        alpha = calcAsymmetry(sLam);
        BinnedEstimatePtr<string> h_lam;
        book(h_lam, 1, 1, 2);
       	h_lam->bin(1).set(alpha.first, alpha.second);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr  _h_Kp,_h_Km,_h_lm,_h_lb;
    const string Ecm = "91.2";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_1995_I382285);


}
