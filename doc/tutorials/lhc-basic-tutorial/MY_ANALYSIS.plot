BEGIN PLOT /MY_ANALYSIS/mass_*
Title=Dilepton invariant mass spectrum
XLabel=$m_{\ell\ell}$ [GeV]
YLabel=$\text{d}\sigma / \text{d} m_{\ell\ell}$ [pb / GeV]
END PLOT

BEGIN PLOT /MY_ANALYSIS/.*bare.*
Title=Dilepton invariant mass spectrum, bare level
END PLOT

BEGIN PLOT /MY_ANALYSIS/.*dressed.*
Title=Dilepton invariant mass spectrum, dressed level
END PLOT

BEGIN PLOT /MY_ANALYSIS/jets_*
Title=Exclusive jet multiplicity
XLabel=$N_\text{jets}$
YLabel=$\text{d}\sigma / \text{d} N_\text{jets}$ [pb]
END PLOT

BEGIN PLOT /MY_ANALYSIS/bjets_*
Title=Exclusive $b$-jet multiplicity
XLabel=$N_b$
YLabel=$\text{d}\sigma / \text{d} N_b$ [pb]
END PLOT

BEGIN PLOT /MY_ANALYSIS/jet1_pT
Title=Leading jet transverse momentum
XLabel=$p_\text{T}$ [GeV]
YLabel=$\text{d}\sigma / \text{d} p_\text{T}$ [pb / GeV]
END PLOT

BEGIN PLOT /MY_ANALYSIS/HT
Title=Scalar jet-$p_\text{T}$ sum
XLabel=$H_\text{T}$ [GeV]
YLabel=$\text{d}\sigma / \text{d} H_\text{T}$ [pb / GeV]
END PLOT

BEGIN PLOT /MY_ANALYSIS/pTmiss
Title=Missing transverse momentum
XLabel=$p_\text{T}^\text{miss}$ [GeV]
YLabel=$\text{d}\sigma / \text{d} p_\text{T}^\text{miss}$ [pb / GeV]
END PLOT
