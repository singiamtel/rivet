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


BEGIN PLOT /EXAMPLE_SIMPLEANALYSIS/count_pass
Title=Signal region event yield
YLabel=Expected event yield
END PLOT
