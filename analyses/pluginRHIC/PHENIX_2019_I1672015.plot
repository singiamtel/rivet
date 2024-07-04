# BEGIN PLOT /PHENIX_2019_I1672015/DiMuon_pT
YMin=1E-5
YMax=10
Title=Drell-Yan production at $\sqrt{s}=200$ GeV
XLabel=$p_T$ (GeV)
YLabel=$1/\sigma d\sigma/dp_T$ (GeV$^{-1}$)
RatioPlot=0

# END PLOT
BEGIN PLOT /PHENIX_2019_I1672015/d01-x01-y01
Title=Drell-Yan $\sqrt{s}=200$ GeV, $ 1.2 < |y| < 2.2$
XLabel=$p_T$ (GeV)
YLabel=$1/(2 \pi p_T) d^2 \sigma/ dy /dp_T$ (pb/GeV$^2$)
LogY=1
LegendTitle=$4.8 < m_{\mu^+\mu^-} < 8.2 $ GeV
YMin=1E-2
YMax=300
RatioPlotYMin=0.01
RatioPlotYMax=4.0
GofType=chi2
END PLOT
BEGIN PLOT /PHENIX_2019_I1672015/d02-x01-y01
Title=PHENIX; Drell-Yan $\sqrt{s}=200$ GeV
YMin=0.5
YMax=500
XLabel=$M$ (GeV)
YLabel=$d^2 \sigma/ dy /dM$ (pb/GeV)
LogY=1
LegendTitle=$ 1.2 < |y| < 2.2$
RatioPlotYMin=0.01
RatioPlotYMax=4.0
END PLOT
BEGIN PLOT /PHENIX_2019_I1672015/DiMuon_y
Title=PHENIX $DY \to \mu^+ \mu^-$, $\sqrt{s} = 200$ GeV,
XLabel=$y$ [GeV]
YLabel=$d \sigma/ dy$
LogY=1
LegendTitle=$4.8 < m_{\mu^+\mu^-} < 8.2 $ GeV
RatioPlotYMin=0.01
RatioPlotYMax=4.0
END PLOT

# ... add more histograms as you need them ...
