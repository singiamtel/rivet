# specific cases where the PDG string can't be parsed
special_modes={}
special_latex={}
# D+ decays
special_modes[411] = { 277 : [r"D+ --> K0S pi+ pi0  nonresonant ",
                              r"D+ --> kappabar()0 pi+, kappabar()0 --> K0S pi0"] }
special_latex[411] = { 277 : r"$D^+\to K^0_S \pi^+\pi^0 (\text{non-res}\ +\ \kappa^0)$"}
# D_s decays
special_modes[431] = { 85 : [r"D_s()+ --> eta e+ nu_e",
                             r"D_s()+ --> eta^'(958) e+ nu_e"],
                       176 : [r"D_s()+ --> Kbar_1(1270)0 K+ , Kbar_1(1270)0 --> K^*(892)- pi+",
                              r"D_s()+ --> Kbar_1(1270)0 K+ , Kbar_1(1270)0 --> Kbar^*(892)0 pi0"],
                       178 : [r"D_s()+ --> Kbar_1(1400)0 K+ , Kbar_1(1400)0 --> K^*(892)- pi+",
                              r"D_s()+ --> Kbar_1(1400)0 K+ , Kbar_1(1400)0 --> Kbar^*(892)0 pi0"],
                       179 : [r"D_s()+  --> eta(1475) pi+ , eta --> K^*(892)0 K0S , K^*()0 --> K- pi+"],
                       183 : [r"D_s()+ --> 2pi+ pi- eta non-eta^'(958)"],
                       199 : [r"D_s()+  --> eta(1475) pi+ , eta --> Kbar^*(892)0 K0S , Kbar^*()0 --> K- pi+"],
                       211 : [r"D_s()+  --> K^*(892)+ K0S,  K^*(892)+  --> K0S pi+"]}
special_latex[431] = { 85  : r"$D_{s}^{+} \to \eta/\eta^{\prime}(958) e^{+} \nu_{e} $",
                       176 : r"$D_{s}^{+} \to \overline{K}_{1}(1270)^{0} (\to \overline{K}^{*}(892) \pi)K^{+}$",
                       178 : r"$D_{s}^{+} \to \overline{K}_{1}(1400)^{0} (\to \overline{K}^{*}(892) \pi)K^{+}$"}
# eta_c decays
special_modes[441] = { 14  : [r"eta_c(1S) --> K+ Kbar0 pi-","eta_c(1S) --> K- K0 pi+",
                              r"eta_c(1S) --> K+ K- pi0","eta_c(1S) --> K0 Kbar0 pi0"],
                       18  : [r"eta_c(1S) --> K^*(892)+ K^*(892)-","eta_c(1S) --> K^*(892)0 Kbar^*(892)0"],
                       19  : [r"eta_c(1S) --> rho+ rho-","eta_c(1S) --> rho0 rho0"],
                       24  : [r"eta_c(1S) --> eta^'(958) pi+ pi-","eta_c(1S) --> eta^'(958) pi0 pi0"],
                       25  : [r"eta_c(1S) --> K0 Kbar0 eta","eta_c(1S) --> K+ K- eta"],
                       85  : [r"eta_c(1S) --> K0 Kbar0 eta^'(958)","eta_c(1S) --> K+ K- eta^'(958)"]}
special_latex[441] = { 14  : r"$\eta_{c}(1S) \to K \overline{K} \pi$",
                       18  : r"$\eta_{c}(1S) \to K^{*}(892) \overline{K}^{*}(892)$",
                       19  : r"$\eta_{c}(1S) \to \rho \rho$",
                       24  : r"$\eta_{c}(1S) \to \eta^{\prime}(958) \pi\pi $",
                       25  : r"$\eta_{c}(1S) \to \eta K \overline{K}$",
                       85  : r"$\eta_{c}(1S) \to \eta^{\prime}(958) K \overline{K}$"}
# J/psi decays
special_modes[443] = { 232 : [r"J/psi(1S) --> K+ K^*(892)- --> Kbar0 K+ pi-",
                              r"J/psi(1S) --> K- K^*(892)+ --> K0 K- pi+"],
                       233 : [r"J/psi(1S) --> K0 Kbar^*(892)0 --> K0 K- pi+",
                              r"J/psi(1S) --> Kbar0 K^*(892)0 --> Kbar0 K+ pi-"]}
special_latex[443] = { 176 : r"$D_{s}^{+} \to \overline{K}_{1}(1270)^{0} (\to \overline{K}^{*}(892) \pi )K^{+} $",
                       178 : r"$D_{s}^{+} \to \overline{K}_{1}(1400)^{0} (\to \overline{K}^{*}(892) \pi )K^{+} $",
                       232 : r"$J/\psi(1S) \to K^+ K^*(892)^- \to \bar{K}^0 K^+ \pi^- (+\text{c.c.})$",
                       233 : r"$J/\psi(1S) \to K^0 \bar{K}^*(892)^0 \to K^0 K^- \pi^+ (+\text{c.c.})$"}
