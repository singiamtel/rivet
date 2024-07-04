# BEGIN PLOT /R209_1982_I168182/DiMuon_pT
YMin=1E-5
YMax=10
Title=Drell-Yan production at $\sqrt{s}=62$ GeV
XLabel=$p_T$ (GeV)
YLabel=$1/\sigma d\sigma/dp_T$ (GeV$^{-1}$)
RatioPlot=0

# END PLOT
BEGIN PLOT /R209_1982_I168182/d0*
Title=Drell-Yan $\sqrt{s}=62$ GeV
END PLOT
BEGIN PLOT /R209_1982_I168182/d01-x01-y01
YMin=5E-5
YMax=5
XMax=20
XLabel=$M$ (GeV)
YLabel=$d\sigma/dM$ (bb/GeV)
END PLOT
BEGIN PLOT /R209_1982_I168182/d01-x01-y02
LegendTitle=$3.5 < M_{\mu^+\mu^-} < 18$ GeV
XLabel=$M$ (GeV)
YLabel=$d\sigma/dM$(nb/GeV)
END PLOT
BEGIN PLOT /R209_1982_I168182/d02-x01-y01
LegendTitle=$5 < M_{\mu^+\mu^-} < 8$ GeV
YMin=1E-4
YMax=5
XLabel=$p_T$ (GeV)
YLabel=$1/(2 p_T) d\sigma/dp_T$ (nb/GeV$^2$)
YMax=1
YMin=1E-4
GofType=chi2
END PLOT

# ... add more histograms as you need them ...
