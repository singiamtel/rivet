#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {


  /// ttbar lepton+jets at 13 TeV
  class CMS_2021_I1901295 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1901295);

    const double mtop_s = 172.5;
    const double mw_s = 80.4;

    void init() {

      const FinalState fs(Cuts::abseta < 6.);

      declare(fs, "FS");

      declare(InvisibleFinalState(), "Invisibles");

      FinalState all_photons(Cuts::abspid == PID::PHOTON);

      PromptFinalState prompt_leptons(Cuts::abseta < 6. && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);

      LeptonFinder dressed_leptons(prompt_leptons, all_photons, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 15*GeV);
      declare(dressed_leptons, "MyLeptons");

      declare(FastJets(fs, JetAlg::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.8), "JetsAK8");

      book(_h["thadpt"], "d159-x01-y01");
      book(_h["tleppt"], "d163-x01-y01");
      book(_h["thardpt"], "d167-x01-y01");
      book(_h["tsoftpt"], "d171-x01-y01");
      book(_h["st"], "d175-x01-y01");
      book(_h["thady"], "d179-x01-y01");
      book(_h["tlepy"], "d183-x01-y01");
      book(_h["dy"], "d187-x01-y01");
      book(_h["ady"], "d191-x01-y01");
      book(_h["ttm"], "d195-x01-y01");
      book(_h["ttpt"], "d199-x01-y01");
      book(_h["tty"], "d203-x01-y01");
      book(_h["ttphi"], "d207-x01-y01");
      book(_h["cts"], "d211-x01-y01");
      book(_h["lpt"], "d317-x01-y01");
      book(_h["njet"], "d321-x01-y01");
      book(_h["ht"], "d325-x01-y01");
      book(_h["mevt"], "d329-x01-y01");

      book(_h["thadpt_norm"], "d161-x01-y01");
      book(_h["tleppt_norm"], "d165-x01-y01");
      book(_h["thardpt_norm"], "d169-x01-y01");
      book(_h["tsoftpt_norm"], "d173-x01-y01");
      book(_h["st_norm"], "d177-x01-y01");
      book(_h["thady_norm"], "d181-x01-y01");
      book(_h["tlepy_norm"], "d185-x01-y01");
      book(_h["dy_norm"], "d189-x01-y01");
      book(_h["ady_norm"], "d193-x01-y01");
      book(_h["ttm_norm"], "d197-x01-y01");
      book(_h["ttpt_norm"], "d201-x01-y01");
      book(_h["tty_norm"], "d205-x01-y01");
      book(_h["ttphi_norm"], "d209-x01-y01");
      book(_h["cts_norm"], "d213-x01-y01");
      book(_h["lpt_norm"], "d319-x01-y01");
      book(_h["njet_norm"], "d323-x01-y01");
      book(_h["ht_norm"], "d327-x01-y01");
      book(_h["mevt_norm"], "d331-x01-y01");

      vector<double> thadptbins = {0.0, 80.0, 160.0, 240.0, 320.0, 400.0, 1500.0};
      book(_b["thadpt_thady"], thadptbins);
      book(_b["thadpt_thady_norm"], thadptbins);
      for (size_t i = 0; i < _b["thadpt_thady"]->numBins(); ++i) {
        book(_b["thadpt_thady"]->bin(i+1), 215 + i, 1, 1);
        book(_b["thadpt_thady_norm"]->bin(i+1), 222 + i, 1, 1);
      }
      vector<double> ttmbins = {250.0, 420.0, 520.0, 620.0, 800.0, 1000.0, 3500.0};
      book(_b["ttm_tty"], ttmbins);
      book(_b["ttm_tty_norm"], ttmbins);
      book(_b["ttm_cts"], ttmbins);
      book(_b["ttm_cts_norm"], ttmbins);
      book(_b["ttm_thadpt"], ttmbins);
      book(_b["ttm_thadpt_norm"], ttmbins);
      for (size_t i = 0; i < _b["ttm_tty"]->numBins(); ++i) {
        book(_b["ttm_tty"]->bin(i+1), 229 + i, 1, 1);
        book(_b["ttm_tty_norm"]->bin(i+1), 236 + i, 1, 1);
        book(_b["ttm_cts"]->bin(i+1), 243 + i, 1, 1);
        book(_b["ttm_cts_norm"]->bin(i+1), 250 + i, 1, 1);
        book(_b["ttm_thadpt"]->bin(i+1), 257 + i, 1, 1);
        book(_b["ttm_thadpt_norm"]->bin(i+1), 264 + i, 1, 1);
      }
      vector<double> ttptbins = {0.0, 50.0, 120.0, 200.0, 300.0, 400.0, 1200.0};
      book(_b["ttpt_thadpt"], ttptbins);
      book(_b["ttpt_thadpt_norm"], ttptbins);
      for (size_t i = 0; i < _b["ttpt_thadpt"]->numBins(); ++i) {
        book(_b["ttpt_thadpt"]->bin(i+1), 271 + i, 1, 1);
        book(_b["ttpt_thadpt_norm"]->bin(i+1), 278 + i, 1, 1);
      }
      vector<double> adybins = {0.0, 0.6, 1.2, 1.8, 3.5};
      book(_b["ady_ttm"], adybins);
      book(_b["ady_ttm_norm"], adybins);
      for (size_t i = 0; i < _b["ady_ttm"]->numBins(); ++i) {
        book(_b["ady_ttm"]->bin(i+1), 285 + i, 1, 1);
        book(_b["ady_ttm_norm"]->bin(i+1), 290 + i, 1, 1);
      }
      vector<double> ttmbins2 = {250.0, 700.0, 1000.0, 3500.0};
      book(_b["ttm_dy"], ttmbins2);
      book(_b["ttm_dy_norm"], ttmbins2);
      for (size_t i = 0; i < _b["ttm_dy"]->numBins(); ++i) {
        book(_b["ttm_dy"]->bin(i+1), 295 + i, 1, 1);
        book(_b["ttm_dy_norm"]->bin(i+1), 299 + i, 1, 1);
      }
      vector<double> topybins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
      book(_b["topy_topbary"], topybins);
      book(_b["topy_topbary_norm"], topybins);
      for (size_t i = 0; i < _b["topy_topbary"]->numBins(); ++i) {
        book(_b["topy_topbary"]->bin(i+1), 303 + i, 1, 1);
        book(_b["topy_topbary_norm"]->bin(i+1), 311 + i, 1, 1);
      }
      vector<double> njetbins = {-0.5, 0.5, 1.5, 2.5, 3.5};
      book(_b["njet_ttm"], njetbins);
      book(_b["njet_ttm_norm"], njetbins);
      book(_b["njet_ttpt"], njetbins);
      book(_b["njet_ttpt_norm"], njetbins);
      book(_b["njet_thadpt"], njetbins);
      book(_b["njet_thadpt_norm"], njetbins);
      for (size_t i = 0; i < _b["njet_ttm"]->numBins(); ++i) {
        book(_b["njet_ttm"]->bin(i+1), 333 + i, 1, 1);
        book(_b["njet_ttm_norm"]->bin(i+1), 338 + i, 1, 1);
        book(_b["njet_ttpt"]->bin(i+1), 343 + i, 1, 1);
        book(_b["njet_ttpt_norm"]->bin(i+1), 348 + i, 1, 1);
        book(_b["njet_thadpt"]->bin(i+1), 353 + i, 1, 1);
        book(_b["njet_thadpt_norm"]->bin(i+1), 358 + i, 1, 1);
      }
    }

    void analyze(const Event &event) {
      DressedLeptons m_leptons;
      Jets m_bjets;
      Jets m_ljets;
      Jets m_alljets;
      Particles m_thadboosted;
      Particles m_tlepboosted;
      Particles m_additionaljets;

      int numvetoleps = 0;
      const DressedLeptons &dressedleptons =
          apply<LeptonFinder>(event, "MyLeptons").dressedLeptons();
      for (const DressedLepton &lep : dressedleptons) {
        if (lep.pt()/GeV > 30. && lep.abseta() < 2.4) {
          m_leptons.push_back(lep);
        } else {
          ++numvetoleps;
        }
      }
      if (m_leptons.size() != 1 || numvetoleps != 0) {
        vetoEvent;
      }

      DressedLepton lepton = m_leptons[0];

      const Particles &invfspars =
          apply<FinalState>(event, "Invisibles").particles(Cuts::abseta < 6.);
      FourMomentum nusum = accumulate(
          invfspars.begin(), invfspars.end(), FourMomentum(0., 0., 0., 0.),
          [&](const FourMomentum &invmom, const Particle &par) {
            return invmom + par.momentum();
          });

      const Jets &allAK4Jets =
          apply<FastJets>(event, "JetsAK4")
              .jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 25. * GeV);

      for (const Jet &jet : allAK4Jets) {
        if (deltaR(lepton, jet) > 0.4) {
          if (jet.bTagged()) {
            m_bjets.push_back(jet);
          } else {
            m_ljets.push_back(jet);
          }

          m_alljets.push_back(jet);
        } else if (jet.bTagged() && lepton.pt()/GeV > 50.) {
          const Particle &undressedlep(lepton.bareLepton());
          bool injet = false;
          for (const Particle &con : jet.particles()) {
            if (con.momentum() == undressedlep.momentum()) {
              injet = true;
              break;
            }
          }
          FourMomentum purejet = jet.momentum();
          if (injet) {
            purejet -= lepton.momentum();
          }
          double checkreco;
          FourMomentum nureco =
              numom(nusum, lepton.momentum(), purejet, checkreco);
          if (checkreco < 0.) {
            continue;
          }
          FourMomentum tlepmom(nureco + purejet + lepton.momentum());
          if (!injet) {
            tlepmom += lepton.momentum();
          }
          if (tlepmom.pt()/GeV > 380. && tlepmom.mass()/GeV > 100. &&
              tlepmom.mass()/GeV < 240.) {
            m_tlepboosted.emplace_back(6, tlepmom);
          }
        }
      }

      const Jets &allAK8Jets =
          apply<FastJets>(event, "JetsAK8")
              .jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 380. * GeV);
      for (const Jet &jet : allAK8Jets) {
        if (deltaR(lepton, jet) < 0.8) {
          continue;
        }
        if (jet.bTagged() && jet.mass()/GeV > 120.) {
          m_thadboosted.emplace_back(6, jet.momentum());
        }
      }

      bool reco = false;
      Particle thad;
      Particle tlep;

      // boosted thad, boosted tlep reconstruction
      if (m_tlepboosted.size() != 0) {
        double bestdm = numeric_limits<double>::max();
        sort(m_tlepboosted.begin(), m_tlepboosted.end(),
             [&mtop = mtop_s](const FourMomentum &A, const FourMomentum &B) {
               return abs(A.mass()/GeV - mtop) < abs(B.mass()/GeV - mtop);
             });
        tlep = m_tlepboosted[0];
        for (const Particle &th : m_thadboosted) {
          if (deltaR(th, tlep) < 1.2) {
            continue;
          }
          double dm = abs(th.mass()/GeV - mtop_s);
          if (dm < bestdm) {
            bestdm = dm;
            thad = Particle(6, th.momentum());
            reco = true;
          }
        }
      }

      if (reco) {
        for (const Jet &jet : m_alljets) {
          if (deltaR(thad, jet) < 1.2) {
            continue;
          }
          if (deltaR(tlep, jet) < 0.6) {
            continue;
          }
          m_additionaljets.emplace_back(1, jet.momentum());
        }
      } else {
        // boosted thad, resolved tlep reconstruction
        double bestcoma = numeric_limits<double>::max();
        Particle tlepcoma;
        Particle thadcoma;
        Particle blepcoma;
        for (const Particle &th : m_thadboosted) {
          for (const Jet &bjl : m_bjets) {
            if (deltaR(th, bjl) < 1.2) {
              continue;
            }
            if (deltaR(th, lepton) < 1.2) {
              continue;
            }
            double checkreco;
            FourMomentum nureco = numom(nusum, lepton.momentum(), bjl, checkreco);
            if (checkreco < 0.) {
              continue;
            }
            Particle tl =
                Particle(6, nureco + lepton.momentum() + bjl.momentum());
            if (tl.mass()/GeV < 100. || tl.mass()/GeV > 240.) {
              continue;
            }
            double coma = pow(th.mass()/GeV - mtop_s, 2) + pow(tl.mass()/GeV - mtop_s, 2);
            if (coma < bestcoma) {
              bestcoma = coma;
              blepcoma = Particle(5, bjl.momentum());
              //!!!To reproduce the distributions in the paper, we have to use
              //!nusum as neutrino momentum here instead of nureco.!!!
              tlepcoma = Particle(6, nusum + lepton.momentum() + bjl.momentum());
              // tlepcoma = Particle(6, tl);
              thadcoma = th;
            }
          }
        }

        double bestcomb = numeric_limits<double>::max();
        Particle tlepcomb;
        Particle thadcomb;
        Particles ttdecay(4);
        if (m_bjets.size() >= 2 && m_ljets.size() >= 2) {
          // resolved thad, resolved tlep reconstruction
          for (const Jet &bjl : m_bjets) {
            double checkreco;
            FourMomentum nureco =
                numom(nusum, lepton.momentum(), bjl.momentum(), checkreco);
            if (checkreco < 0.) {
              continue;
            }
            FourMomentum tl(lepton.momentum() + nureco + bjl.momentum());
            if (tl.mass()/GeV < 100. || tl.mass()/GeV > 240.) {
              continue;
            }

            for (size_t a = 0; a < m_ljets.size(); ++a) {
              const Jet &lja = m_ljets[a];
              for (size_t b = 0; b < a; ++b) {
                const Jet &ljb = m_ljets[b];
                FourMomentum wh(lja.momentum() + ljb.momentum());
                for (const Jet &bjh : m_bjets) {
                  if (&bjh == &bjl) {
                    continue;
                  }
                  FourMomentum th(wh + bjh.momentum());
                  if (th.mass()/GeV < 100. || th.mass()/GeV > 240.) {
                    continue;
                  }

                  double comb = pow(wh.mass()/GeV - mw_s, 2) +
                                pow(th.mass()/GeV - mtop_s, 2) +
                                pow(tl.mass()/GeV - mtop_s, 2);
                  if (comb < bestcomb) {
                    bestcomb = comb;
                    thadcomb = Particle(6, th);
                    //!!!To reproduce the distributions in the paper, we have to
                    //!use nusum as neutrino momentum here instead of nureco.!!!
                    tlepcomb =
                        Particle(6, lepton.momentum() + bjl.momentum() + nusum);
                    // tlepcomb = Particle(6, tl);
                    ttdecay[0] = Particle(5, bjh);
                    ttdecay[1] = Particle(1, lja);
                    ttdecay[2] = Particle(1, ljb);
                    ttdecay[3] = Particle(5, bjl);
                  }
                }
              }
            }
          }
        }

        if (bestcoma != numeric_limits<double>::max() &&
            bestcomb != numeric_limits<double>::max()) {
          if (abs(thadcoma.mass()/GeV - mtop_s) < abs(thadcomb.mass()/GeV - mtop_s)) {
            bestcomb = numeric_limits<double>::max();
          } else {
            bestcoma = numeric_limits<double>::max();
          }
        }

        if (bestcoma != numeric_limits<double>::max()) {
          reco = true;
          thad = thadcoma;
          tlep = tlepcoma;
          for (const Jet &jet : m_alljets) {
            if (deltaR(blepcoma, jet) < 0.01) {
              continue;
            }
            if (deltaR(thadcoma, jet) < 1.2) {
              continue;
            }
            m_additionaljets.emplace_back(1, jet.momentum());
          }
        } else if (bestcomb != numeric_limits<double>::max()) {
          reco = true;
          thad = thadcomb;
          tlep = tlepcomb;
          for (const Jet &jet : m_alljets) {
            if (find_if(ttdecay.begin(), ttdecay.end(), [&](const Particle &par) {
                  return deltaR(jet, par) < 0.01;
                }) != ttdecay.end()) {
              continue;
            }
            m_additionaljets.emplace_back(1, jet.momentum());
          }
        }
      }

      if (!reco) {
        vetoEvent;
      }
      FourMomentum tt(thad.momentum() + tlep.momentum());
      FourMomentum top = lepton.pid() > 0 ? thad : tlep;
      FourMomentum topbar = lepton.pid() < 0 ? thad : tlep;

      const double ht = sum(m_additionaljets, Kin::pT, 0.)/GeV;
      FourMomentum mevt =  std::accumulate(m_additionaljets.begin(), m_additionaljets.end(), tt,
                     [](const FourMomentum &mevt, const Particle &jet) {
                       return mevt + jet;
                     });

      LorentzTransform boostcms;
      boostcms.setBetaVec(-1. * tt.betaVec());
      FourMomentum topcms = boostcms(top);
      double cts = topcms.vector3().dot(tt.vector3()) / (topcms.p() * tt.p());
      double day = abs(top.rapidity()) - abs(topbar.rapidity());
      double ady = abs(top.rapidity() - topbar.rapidity());
      double ttphi = abs(deltaPhi(top, topbar)) * 180. / M_PI;

      dualfill("thadpt", thad.pt()/GeV);
      dualfill("tleppt", tlep.pt()/GeV);
      dualfill("thardpt", max(thad.pt()/GeV, tlep.pt()/GeV));
      dualfill("tsoftpt", min(thad.pt()/GeV, tlep.pt()/GeV));
      dualfill("st", thad.pt()/GeV + tlep.pt()/GeV);
      dualfill("thady", abs(thad.rapidity()));
      dualfill("tlepy", abs(tlep.rapidity()));
      dualfill("dy", day);
      dualfill("ady", ady);
      dualfill("ttm", tt.mass()/GeV);
      dualfill("ttpt", tt.pt()/GeV);
      dualfill("tty", abs(tt.rapidity()));
      dualfill("ttphi", ttphi);
      dualfill("cts", cts);
      dualfill("lpt", lepton.pt()/GeV);
      dualfill("njet", min(m_additionaljets.size(), 6u) + 0.5);
      dualfill("ht", ht/GeV);
      dualfill("mevt", mevt.mass()/GeV);

      dualfill("thadpt_thady", thad.pt()/GeV, abs(thad.rapidity()));
      dualfill("ttm_tty", tt.mass()/GeV, abs(tt.rapidity()));
      dualfill("ttm_cts", tt.mass()/GeV, cts);
      dualfill("ttm_thadpt", tt.mass()/GeV, thad.pt()/GeV);
      dualfill("ttpt_thadpt", tt.pt()/GeV, thad.pt()/GeV);
      dualfill("ady_ttm", ady, tt.mass()/GeV);
      dualfill("ttm_dy", tt.mass()/GeV, day);
      dualfill("topy_topbary", abs(top.rapidity()), abs(topbar.rapidity()));
      dualfill("njet_ttm", min(m_additionaljets.size(), 3u), tt.mass()/GeV);
      dualfill("njet_ttpt", min(m_additionaljets.size(), 3u), tt.pt()/GeV);
      dualfill("njet_thadpt", min(m_additionaljets.size(), 3u), thad.pt()/GeV);
    }

    void finalize() {

      const double sf = crossSection()/picobarn/sumOfWeights();

      for (auto& item : _h) {
        if (item.first.find("_norm") != string::npos) {
          normalize(item.second, 1.0, false);
        }
        else  scale(item.second, sf);
      }

      for (auto& item : _b) {
        if (item.first.find("_norm") != string::npos) {
          normalizeGroup(item.second, 1.0, false);
        }
        else {
          scale(item.second, sf);
        }
        divByGroupWidth(item.second);
      }
    }

    map<string,Histo1DPtr> _h;
    map<string,Histo1DGroupPtr> _b;

    void dualfill(const string& tag, const double value) {
      _h[tag]->fill(value);  _h[tag + "_norm"]->fill(value);
    }

    void dualfill(const string& tag, const double val1, const double val2) {
      _b[tag]->fill(val1, val2);  _b[tag + "_norm"]->fill(val1, val2);
    }

    FourMomentum numom(const FourMomentum &met, const FourMomentum &l,
                              const FourMomentum &b, double &mtres) {
      mtres = -1.;
      FourMomentum va;
      FourMomentum vb;

      double tltn = l.px() * met.px() + l.py() * met.py();
      double tntn = met.px() * met.px() + met.py() * met.py();

      double C = pow(l.pz() / l.E(), 2) - 1.;
      double B =
          pow(mw_s / l.E(), 2) * l.pz() + 2. * l.pz() / l.E() / l.E() * tltn;
      double A = -1. * tntn + pow(mw_s * mw_s / 2. / l.E(), 2) +
                 pow(mw_s / l.E(), 2) * tltn + pow(tltn / l.E(), 2);

      double D = pow(B / C / 2., 2) - A / C;
      if (D >= 0.) {
        va.setXYZM(met.px(), met.py(), -0.5 * B / C + sqrt(D), 0.);
        vb.setXYZM(met.px(), met.py(), -0.5 * B / C - sqrt(D), 0.);
        double ma = (va + l + b).mass()/GeV;
        double mb = (vb + l + b).mass()/GeV;

        if (abs(ma - mtop_s) < abs(mb - mtop_s)) {
          mtres = ma;
          return va;
        } else {
          mtres = mb;
          return vb;
        }
      } else {
        double As = 0.25 * (pow(mw_s * mw_s * l.pz() / l.E() / l.E(), 2) / C -
                            pow(mw_s * mw_s / l.E(), 2));
        double Bsx =
            (pow(mw_s * l.pz() / l.E() / l.E(), 2) / C - pow(mw_s / l.E(), 2)) *
            l.px() * met.px();
        double Bsy =
            (pow(mw_s * l.pz() / l.E() / l.E(), 2) / C - pow(mw_s / l.E(), 2)) *
            l.py() * met.py();
        double Csxx = (pow(l.pz() / l.E() / l.E(), 2) / C - 1. / l.E() / l.E()) *
                          l.px() * met.px() * l.px() * met.px() +
                      met.px() * met.px();
        double Csxy = (pow(l.pz() / l.E() / l.E(), 2) / C - 1. / l.E() / l.E()) *
                      2. * l.px() * met.px() * l.py() * met.py();
        double Csyy = (pow(l.pz() / l.E() / l.E(), 2) / C - 1. / l.E() / l.E()) *
                          l.py() * met.py() * l.py() * met.py() +
                      met.py() * met.py();

        As /= C;
        Bsx /= C;
        Bsy /= C;
        Csxx /= C;
        Csxy /= C;
        Csyy /= C;

        double x = 1.;
        double y = 1.;

        double U =
            As + x * Bsx + y * Bsy + x * x * Csxx + x * y * Csxy + y * y * Csyy;

        double step = 0.1;

        while (true) {
          double dx = Bsx + 2. * Csxx * x + Csxy * y;
          double dy = Bsy + 2. * Csyy * y + Csxy * x;

          x += step * dx / sqrt(dx * dx + dy * dy);
          y += step * dy / sqrt(dx * dx + dy * dy);

          double nU = (As + x * Bsx + y * Bsy + x * x * Csxx + x * y * Csxy +
                       y * y * Csyy);

          if (nU * U < 0.) {
            step *= -0.5;
          }
          U = nU;
          if (abs(U) < 0.01) {
            break;
          }
        }
        double pz = -0.5 / C *
                    (pow(mw_s / l.E(), 2) * l.pz() +
                     2. * l.pz() / l.E() / l.E() *
                         (l.px() * met.px() * x + l.py() * met.py() * y));
        va.setXYZM(met.px() * x, met.py() * y, pz, 0.);
        mtres = (va + l + b).mass()/GeV;
        return va;
      }

      return va;
    }

  };

  RIVET_DECLARE_PLUGIN(CMS_2021_I1901295);

} // namespace Rivet
