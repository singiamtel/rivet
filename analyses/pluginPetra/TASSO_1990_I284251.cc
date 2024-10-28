// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief K0 and K*+ spectra
  class TASSO_1990_I284251 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1990_I284251);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Sphericity(cfs), "Sphericity");


      // Book histograms
      _sqs = 10; size_t i = 0;
      for (double eVal : allowedEnergies()) {

        if (isCompatibleWithSqrtS(eVal))  _sqs = i;

        size_t _ih, _iy;
        if (i==0) { _ih=1; _iy=0; }
        else if (i==1) { _ih=2; _iy=0; }
        else if (i==2) { _ih=0; _iy=3; }
        else if (i==3) { _ih=0; _iy=2; }
        else           { _ih=0; _iy=1; }

        if (_ih==0) {
          book(_h_K0_x[i], 1, 1, _iy);
          if (_iy != 3) {
            book(_p_K0_S_1[i], 5, 1, 2*_iy-1);
            book(_p_K0_S_2[i],"TMP/p_K0_S_2_"+to_string(i),refData(5, 1, 2*_iy));
          }
          book(_h_Kstar_x[i-2], 8, 1, _iy);
          if (_iy == 2) {
            book(_p_Kstar_S_1, 10, 1, 1);
            book(_p_Kstar_S_2, "TMP/p_Kstar_S_2", refData(10,1,2));
          }
        }
        else {
          book(_h_K0_x[i],   _ih+1, 1, 1);
          book(_p_K0_S_1[i], _ih+5, 1, 1);
          book(_p_K0_S_2[i], "TMP/p_K0_S_2"+to_string(i), refData(_ih+5,1,2));
        }
        ++i;
      }
      if (_sqs == 10 && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
      book(_n_K0,    4, 1, 1);
      book(_n_Kstar, 9, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
				   beams.second.p3().mod() ) / 2.0;
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      unsigned int nK0(0),nKstar(0);
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==323 or Cuts::pid==130 or Cuts::pid==310)) {
        double xE = p.E()/meanBeamMom;
        double modp = p.p3().mod();
        double beta = modp/p.E();
        if (p.abspid()==323) {
          if (_sqs>1)  _h_Kstar_x[_sqs-2]->fill(xE,1./beta);
          ++nKstar;
        }
        else {
          _h_K0_x[_sqs]->fill(xE,1./beta);
          ++nK0;
        }
      }
      _n_K0   ->fill(_dedges[_sqs],nK0   );
      _n_Kstar->fill(_dedges[_sqs],nKstar);
      const double sphere = sphericity.sphericity();
      if (_sqs != 2) {
        _p_K0_S_1[_sqs]->fill(sphere,nK0);
        _p_K0_S_2[_sqs]->fill(sphere,cfs.particles().size());
      }
      if (_sqs == 3) {
        _p_Kstar_S_1->fill(sphere,nKstar);
        _p_Kstar_S_2->fill(sphere,cfs.particles().size());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = sqr(sqrtS())*crossSection()/sumOfWeights();
      scale(_h_K0_x,   sf/microbarn);
      scale(_h_Kstar_x, sf/nanobarn);

      size_t i=0;
      for (double eVal : allowedEnergies()) {

        if (isCompatibleWithSqrtS( eVal*GeV ))  _sqs = i;

        size_t _ih, _iy;
        if (i==0) { _ih=1; _iy=0; }
        else if (i==1) { _ih=2; _iy=0; }
        else if (i==2) { _ih=0; _iy=3; }
        else if (i==3) { _ih=0; _iy=2; }
        else           { _ih=0; _iy=1; }

        if (_iy != 3) {
          Estimate1DPtr temp;
          if (_ih==0)  book(temp,5,1,2*_iy);
          else         book(temp,_ih+5,1,2);
          divide(_p_K0_S_1[i],_p_K0_S_2[i],temp);
        }
        if (i == 3) {
          Estimate1DPtr temp;
          book(temp,10,1,2);
          divide(_p_Kstar_S_1,_p_Kstar_S_2,temp);
        }
        ++i;
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_K0_x[5], _h_Kstar_x[3];

    Profile1DPtr _p_K0_S_1[5], _p_K0_S_2[5], _p_Kstar_S_1, _p_Kstar_S_2;

    BinnedProfilePtr<string> _n_K0,_n_Kstar;

    vector<string> _dedges{"14.8"s, "21.5"s, "34.5"s, "35.0"s, "42.6"s};

    size_t _sqs;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TASSO_1990_I284251);


}
