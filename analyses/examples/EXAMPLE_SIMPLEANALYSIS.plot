BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/num_.*
XLabel=$N$
YLabel=Efficiency
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/num_elec
Title=Electron multiplicity
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/num_muon
Title=Muon multiplicity
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/num_tau
Title=Hadronic tau multiplicity
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/num_photon
Title=Photon multiplicity
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/num_bjet
Title=$b$-jet multiplicity
END PLOT


BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/pt_.*
XLabel=$p_\mathrm{T}$ [GeV]
YLabel=Acceptance
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/pt_bjet1
Title=Transverse momentum of leading $b$-jet
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/pt_lep1
Title=Transverse momentum of leading charged lepton
END PLOT


BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/met.*
YLabel=Acceptance
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/met
Title=Missing transverse momentum
XLabel=$p_\mathrm{T}^\mathrm{miss}$ [GeV]
END PLOT

BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/met_signf
Title=Significance of missing transverse momentum
XLabel=MET significance $p_\mathrm{T}^\mathrm{miss}/\sigma_\mathrm{MET}$
END PLOT


BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/count_pass
Title=Signal region event yield
YLabel=Expected event yield
END PLOT
