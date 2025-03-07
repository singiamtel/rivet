## Plotting

Rivet provides a `rivet-mkhtml` script for comparing and plotting data files.
This command produces nice comparison plots of publication quality from the YODA files.

The high-level program `rivet-mkhtml` will automatically create a plot
webpage from the given YODA files. It searches for reference data automatically
and uses the other two commands internally. Example:

```
rivet-mkhtml withUE.yoda:'Title=With UE' withoutUE.yoda:'LineColor=blue'
```

Run `rivet-mkhtml --help` to find out about all features and options.

Plotting options will be taken from `*.plot` files installed in Rivet's share directory.
These contain plotting instructions as documented
[here](https://gitlab.com/hepcedar/yoda/-/blob/release-2-1-x/doc/PlotConfig.md).
Such files can also be written for any plugin analysis and will be found if they are in the `RIVET_ANALYSIS_PATH`.

YODA's plotting API is discussed [here](https://gitlab.com/hepcedar/yoda/-/blob/release-2-1-x/doc/PlotAPI.md).


### Changing the reference data style

The style of the reference data can be influenced using a dummy `REF` argument
on the command line, e.g. `REF:LineColor=red`. Additional tags can be supplied
in order to superimpose the reference curve with different properties, e.g.
```
REF:ErrorBandColor=green yellow REF2:ErrorPatter=^stat:ErrorBandColor=darkgreen"
```

### Under- and overflows

Under- and overflow bins aren't plotted by default because their bin width is infinite.
A dummy bin width can be supplied when turning histograms into estimates, e.g. using
```cpp
Estimat1DPtr est = hist->mkEstimate(/*path =*/ hist->path(),
                                    /*source =*/ "stats",
                                    /*divbyvol =*/ true,
                                    /*overflowWidth =*/ 1.0);
```
which will avoid a division by infinity in the overflows.
The `IncludeOverflows=<0|1>` toggle can be used to get these bins visualised:
```
rivet-mkhtml input.yoda.gz PLOT:IncludeOverflows=1
```

