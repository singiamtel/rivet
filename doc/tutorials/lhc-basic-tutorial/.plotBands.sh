
# define custom colours
DBLUE="darkblue"
MBLUE="blue"
LBLUE="lightblue"
# some useful style options for bands from
# https://gitlab.com/hepcedar/yoda/-/blob/release-2-1-x/doc/PlotConfig.md?ref_type=heads#band-styles-append-to-file
MEPS_STYLE="ErrorBandColor=$LBLUE:LineColor=$LBLUE"
ME_STYLE="ErrorBandColor=$MBLUE:LineColor=$MBLUE"
STAT_STYLE="ErrorBandStyle=/:ErrorBandColor=$DBLUE:LineColor=$DBLUE"
######################################################################
# Complete regular expressions here to select the correct muliweights
######################################################################
MEPS_VARS="^MU.*PDF.*"
ME_VARS="ME_ONLY_MU.*PDF.*"
######################################################################
MEPS_BAND="BandComponentEnv=$MEPS_VARS"
ME_BAND="BandComponentEnv=$ME_VARS"
# use "Name" tag declare multiple distinct curves using the same file,
# apply style options, set a title, select the multiweights to be
# combined into bands
rivet-mkhtml -o plots_with_bands3 \
Rivet.yoda":$MEPS_STYLE:$MEPS_BAND:Title=ME+PS scales \$\\oplus\$ stats" \
Rivet.yoda":$ME_STYLE:$ME_BAND:Title=ME scales \$\\oplus\$ stats" \
Rivet.yoda":$STAT_STYLE:Title=stats only"

