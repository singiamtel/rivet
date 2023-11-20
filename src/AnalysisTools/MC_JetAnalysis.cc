// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  MC_JetAnalysis::MC_JetAnalysis(const string& name,
                                 size_t njet,
                                 const string& jetpro_name,
                                 double jetptcut)
    : Analysis(name), _njet(njet), _jetpro_name(jetpro_name), _jetptcut(jetptcut),
      _h_pT_jet(njet),
      _h_eta_jet(njet), _h_eta_jet_plus(njet), _h_eta_jet_minus(njet),
      _h_rap_jet(njet), _h_rap_jet_plus(njet), _h_rap_jet_minus(njet),
      _h_mass_jet(njet), tmpeta(njet), tmprap(njet)
  {    }


  // Book histograms
  void MC_JetAnalysis::init() {
    const double sqrts = sqrtS() ? sqrtS() : 14*TeV;

    // Get an optional rebinning factor from the analysis options
    const int rebin = getOption<int>("REBIN", 1);

    for (size_t i = 0; i < _njet; ++i) {
      const string pTname = "jet_pT_" + to_str(i+1);
      const double pTmax = 1.0/(double(i)+2.0) * sqrts/GeV/2.0;
      const int nbins_pT = 100/(i+1)/rebin;
      if (pTmax > 10) { // Protection aginst logspace exception, needed for LEP
        book(_h_pT_jet[i] ,pTname, logspace(nbins_pT, 10.0, pTmax));
      }

      const string massname = "jet_mass_" + to_str(i+1);
      const double mmax = 100.0;
      const int nbins_m = 100/(i+1)/rebin;
      book(_h_mass_jet[i] ,massname, logspace(nbins_m, 1.0, mmax));

      const string etaname = "jet_eta_" + to_str(i+1);
      book(_h_eta_jet[i] ,etaname, (i > 1 ? 25 : 50)/rebin, -5.0, 5.0);
      book(_h_eta_jet_plus[i], "_" + etaname + "_plus", i > 1 ? 15 : 25, 0, 5);
      book(_h_eta_jet_minus[i], "_" + etaname + "_minus", i > 1 ? 15 : 25, 0, 5);

      const string rapname = "jet_y_" + to_str(i+1);
      book(_h_rap_jet[i] ,rapname, (i > 1 ? 25 : 50)/rebin, -5.0, 5.0);
      book(_h_rap_jet_plus[i], "_" + rapname + "_plus", i > 1 ? 15 : 25, 0, 5);
      book(_h_rap_jet_minus[i], "_" + rapname + "_minus", i > 1 ? 15 : 25, 0, 5);

      book(tmpeta[i], "jet_eta_pmratio_" + to_str(i+1));
      book(tmprap[i], "jet_y_pmratio_" + to_str(i+1));

      for (size_t j = i+1; j < min(size_t(3), _njet); ++j) {
        const std::pair<size_t, size_t> ij = std::make_pair(i, j);
        const string ijstr = to_str(i+1) + to_str(j+1);

        string detaname = "jets_deta_" + ijstr;
        book(_h_deta_jets[ij], detaname, 25/rebin, -5.0, 5.0);

        string dphiname = "jets_dphi_" + ijstr;
        book(_h_dphi_jets[ij], dphiname, 25/rebin, 0.0, M_PI);

        string dRname = "jets_dR_" + ijstr;
        book(_h_dR_jets[ij], dRname, 25/rebin, 0.0, 5.0);
      }
    }

    vector<int> discbins;
    vector<std::string> ratiobins;
    for (size_t i = 0; i < _njet+3; ++i) {
      discbins.push_back(i);
      if (i) {
        const string label = std::to_string(i) + "/" + std::to_string(i-1);
        ratiobins.push_back(label);
      }
    }
    book(_h_jet_multi_exclusive ,"jet_multi_exclusive", discbins);
    book(_h_jet_multi_inclusive ,"jet_multi_inclusive", discbins);
    book(_h_jet_multi_ratio, "jet_multi_ratio", ratiobins);
    book(_h_jet_HT ,"jet_HT", logspace(50/rebin, _jetptcut, sqrts/GeV/2.0));
    book(_h_mjj_jets, "jets_mjj", 40/rebin, 0.0, sqrts/GeV/2.0);
  }


  // Do the analysis
  void MC_JetAnalysis::analyze(const Event & e) {

    const Jets& jets = apply<FastJets>(e, _jetpro_name).jetsByPt(Cuts::pT > _jetptcut);

    for (size_t i = 0; i < _njet; ++i) {
      if (jets.size() < i+1) continue;
      _h_pT_jet[i]->fill(jets[i].pT()/GeV);
      // Check for numerical precision issues with jet masses
      double m2_i = jets[i].mass2();
      if (m2_i < 0) {
        if (m2_i < -1e-4) {
          MSG_WARNING("Jet mass2 is negative: " << m2_i << " GeV^2. "
                      << "Truncating to 0.0, assuming numerical precision is to blame.");
        }
        m2_i = 0.0;
      }

      // Jet mass
      _h_mass_jet[i]->fill(sqrt(m2_i)/GeV);

      // Jet eta
      const double eta_i = jets[i].eta();
      _h_eta_jet[i]->fill(eta_i);
      (eta_i > 0.0 ? _h_eta_jet_plus : _h_eta_jet_minus)[i]->fill(fabs(eta_i));

      // Jet rapidity
      const double rap_i = jets[i].rapidity();
      _h_rap_jet[i]->fill(rap_i);
      (rap_i > 0.0 ? _h_rap_jet_plus : _h_rap_jet_minus)[i]->fill(fabs(rap_i));

      // Inter-jet properties
      for (size_t j = i+1; j < min(size_t(3),_njet); ++j) {
        if (jets.size() < j+1) continue;
        std::pair<size_t, size_t> ij = std::make_pair(i, j);
        double deta = jets[i].eta()-jets[j].eta();
        double dphi = deltaPhi(jets[i].momentum(),jets[j].momentum());
        double dR = deltaR(jets[i].momentum(), jets[j].momentum());
        _h_deta_jets[ij]->fill(deta);
        _h_dphi_jets[ij]->fill(dphi);
        _h_dR_jets[ij]->fill(dR);
      }
    }

    // Multiplicities
    _h_jet_multi_exclusive->fill(jets.size());
    for (size_t i = 0; i < _njet+2; ++i) {
      if (jets.size() >= i) {
        _h_jet_multi_inclusive->fill(i);
      }
    }

    // HT
    double HT = 0.0;
    for (const Jet& jet : jets) {
      HT += jet.pT();
    }
    _h_jet_HT->fill(HT);

    // mjj
    if (jets.size() > 1) {
      double mjj = (jets[0].momentum() + jets[1].momentum()).mass();
      _h_mjj_jets->fill(mjj);
    }
  }


  // Finalize
  void MC_JetAnalysis::finalize() {
    const double scaling = crossSection()/sumOfWeights();
    for (size_t i = 0; i < _njet; ++i) {
      if (_h_pT_jet[i])  scale(_h_pT_jet[i], scaling);
      scale(_h_mass_jet[i], scaling);
      scale(_h_eta_jet[i], scaling);
      scale(_h_rap_jet[i], scaling);

      // Create eta/rapidity ratio plots
      divide(_h_eta_jet_plus[i], _h_eta_jet_minus[i], tmpeta[i]);
      divide(_h_rap_jet_plus[i], _h_rap_jet_minus[i], tmprap[i]);
    }

    // Scale the d{eta,phi,R} histograms
    scale(_h_deta_jets, scaling);
    scale(_h_dphi_jets, scaling);
    scale(_h_dR_jets, scaling);

    // Fill inclusive jet multi ratio
    for (size_t i = 1; i < _h_jet_multi_inclusive->numBins(); ++i) {
      const string label = std::to_string(i) + "/" + std::to_string(i-1);
      if (_h_jet_multi_inclusive->bin(i).sumW() > 0.0) {
        const double ratio = _h_jet_multi_inclusive->bin(i+1).sumW()/_h_jet_multi_inclusive->bin(i).sumW();
        const double relerr_i = _h_jet_multi_inclusive->bin(i).relErrW();
        const double relerr_j = _h_jet_multi_inclusive->bin(i+1).relErrW();
        const double err = ratio * (relerr_i + relerr_j);
        _h_jet_multi_ratio->binAt(label).set(ratio, {-err,err});
      }
    }

    scale(_h_jet_multi_exclusive, scaling);
    scale(_h_jet_multi_inclusive, scaling);
    scale(_h_jet_HT, scaling);
    scale(_h_mjj_jets, scaling);
  }


}
