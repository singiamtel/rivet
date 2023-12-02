
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief inclusive W/Z cross sections at 7 TeV
  class ATLAS_2016_I1502620 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1502620);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the option system
      _mode = 0;
      _runZ = true;
      _runW = true;
      if ( getOption("LMODE") == "EL" ||
           getOption("LMODE") == "ZEL" ||
           getOption("LMODE") == "WEL" )
        _mode = 1;
      if ( getOption("LMODE") == "MU" ||
           getOption("LMODE") == "ZMU" ||
           getOption("LMODE") == "WMU" )
        _mode = 2;
      if ( getOption("LMODE") == "Z" ||
           getOption("LMODE") == "ZEL" ||
           getOption("LMODE") == "ZMU" )
        _runW = false;
      if ( getOption("LMODE") == "W" ||
           getOption("LMODE") == "WEL" ||
           getOption("LMODE") == "WMU" )
        _runZ = false;


      // Initialise and register projections
      declare("MET", MissingMomentum());
      LeptonFinder ef(0.1, Cuts::pT > 25*GeV && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");
      LeptonFinder mf(0.1, Cuts::pT > 25*GeV && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      DileptonFinder zfe(91.2*GeV, 0.1, Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON, Cuts::massIn(46.0*GeV, 150*GeV));
      declare(zfe, "Zee");
      DileptonFinder zfm(91.2*GeV, 0.1, Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON, Cuts::massIn(46.0*GeV, 150*GeV));
      declare(zfm, "Zmm");


      // Book histograms
      if (_runW) {
        book(_h_Wp_eta,  9, 1, 1);
        book(_h_Wm_eta, 10, 1, 1);
        book(_h_W_asym, 35, 1, 1);
      }
      if (_runZ) {
        book(_h_Zcenlow_y_dressed,  11, 1, 1);
        book(_h_Zcenpeak_y_dressed, 12, 1, 1);
        book(_h_Zcenhigh_y_dressed, 13, 1, 1);
        book(_h_Zfwdpeak_y_dressed, 14, 1, 1);
        book(_h_Zfwdhigh_y_dressed, 15, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // W reco, starting with MET
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();

      // Identify the closest-matching l+MET to m == mW
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const Particles es_mtfilt = select(es, [&](const Particle& e){ return mT(e, pmiss) > 40*GeV; });
      const int iefound = closestMatchIndex(es_mtfilt, pmiss, Kin::mass, 80.4*GeV);
      const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
      const Particles mus_mtfilt = select(mus, [&](const Particle& m){ return mT(m, pmiss) > 40*GeV; });
      const int imfound = closestMatchIndex(mus_mtfilt, pmiss, Kin::mass, 80.4*GeV);

      // Require only one W candidate
      if (pmiss.pT() < 25*GeV && (int(iefound >= 0) + int(imfound >= 0)) == 1 && _runW) {
        Particle lep;
        if (_mode != 2 && iefound >= 0) lep = es_mtfilt[iefound];
        else if (_mode != 1 && imfound >= 0) lep = mus_mtfilt[imfound];

        if (lep.charge3() == 3) _h_Wp_eta->fill(lep.abseta());
        else if (lep.charge3() == -3) _h_Wm_eta->fill(lep.abseta());
      }


      // Now the Z stuff.
      const DileptonFinder& zee = apply<DileptonFinder>(event, "Zee");
      const DileptonFinder& zmm = apply<DileptonFinder>(event, "Zmm");
      if ((zee.bosons().size() + zmm.bosons().size()) == 1 && _runZ) {
        Particle zboson;
        if (_mode != 2 && zee.bosons().size() == 1) zboson = zee.boson();
        else if (_mode !=1 && zmm.bosons().size() == 1 ) zboson = zmm.boson();

        const Particles& leptons = zboson.constituents();
        if (leptons.size() > 1) {
          const double zrap  = zboson.absrap();
          const double zmass = zboson.mass();
          const double eta1 = leptons[0].abseta();
          const double eta2 = leptons[1].abseta();

          // Separation into central/forward and three mass bins
          if (eta1 < 2.5 && eta2 < 2.5) {
            if (zmass < 66.0*GeV)        _h_Zcenlow_y_dressed->fill(zrap);
            else if (zmass < 116.0*GeV)  _h_Zcenpeak_y_dressed->fill(zrap);
            else                         _h_Zcenhigh_y_dressed->fill(zrap);
          }
          else if ((eta1 < 2.5 && 2.5 < eta2 && eta2 < 4.9) || (eta2 < 2.5 && 2.5 < eta1 && eta1 < 4.9)) {
            if (zmass > 66.0*GeV) {
              if (zmass < 116.0*GeV)  _h_Zfwdpeak_y_dressed->fill(zrap);
              else                    _h_Zfwdhigh_y_dressed->fill(zrap);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta)
      if (_runW) asymm(_h_Wp_eta, _h_Wm_eta, _h_W_asym);

      ///  Normalise, scale and otherwise manipulate histograms here
      double lfac = 1.0;
      // If we running on both electrons and muons, divide by two -> xsec for one flavour
      if (_mode == 0) lfac = 0.5;
      const double sf = lfac * 0.5 * crossSection() /picobarn / sumOfWeights(); //< 0.5 accounts for rapidity bin width

      if (_runW){
        scale(_h_Wp_eta, sf);
        scale(_h_Wm_eta, sf);
      }

      if (_runZ){
        scale(_h_Zcenlow_y_dressed, sf);
        scale(_h_Zcenpeak_y_dressed, sf);
        scale(_h_Zcenhigh_y_dressed, sf);
        scale(_h_Zfwdpeak_y_dressed, sf);
        scale(_h_Zfwdhigh_y_dressed, sf);
      }
    }

    /// @}


  protected:

    size_t _mode;
    bool _runZ, _runW;


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_Wp_eta, _h_Wm_eta;
    Estimate1DPtr _h_W_asym;

    Histo1DPtr _h_Zcenlow_y_dressed;
    Histo1DPtr _h_Zcenpeak_y_dressed;
    Histo1DPtr _h_Zcenhigh_y_dressed;
    Histo1DPtr _h_Zfwdpeak_y_dressed;
    Histo1DPtr _h_Zfwdhigh_y_dressed;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1502620);

}
