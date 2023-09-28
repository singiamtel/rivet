
# Bootstrapping event fills

After processing an event, the `AnalysisHandler` can be queried for the fill indices
and fill fractions in order to allow for correlating event fills across various
histograms from multiple analsyses, e.g.

```cpp
print(ah.fillOutcomes())
print(ah.fillFractions())
```

If a file name is supplied to the `AnalysisHandler` via the `setBootstrapFilename` method
this information will be dumped into a corresponding ASCII file.

**Note that this file is updated for every event and will grow in size quickly!**

The format is as follows:

* two header lines at the top of the file (preceded by a `#`): one specifying the weight names
  and one the AO path name,immediately followed by a comma and the total number of bins of the AO
  (including under/overflows and masked bins)
* three lines per event:
  * one starting with a `W` listing the (multi-)weights
  * one starting with an `O` listing the binary fill outcome
  * one starting with an `F` listing the fill fractions

An example using a custom routine that only books a single 1D histogram with four (visible) bins
could look like this:

```
# MEWeight NTrials Weight WeightNormalisation
# /MY_TEST/testHisto,6
W 4.788e-12 5 0.0382951 47.7771
O 0 1 0 1 0 0
F 0 1 0 1 0 0
W 3.26928e-13 16 19.8212 47.7771
O 0 1 0 0 0 0
F 0 2 0 0 0 0
W 6.4679e-11 11 20.3138 47.7771
O 0 1 0 0 1 0
F 0 1 0 0 1 0
```

Note that in the second event, the first visible bin was filled twice.

