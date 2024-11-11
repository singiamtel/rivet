# Merging Histograms 101

There are various reasons why one might want to merge together `YODA` files.
For example, you may want to combine predictions for several sub-processes,
or outputs from multiple jobs for the same process.
Unfortunately, the devil is in the detail, and simply adding the files
is often not enough to get it right.
Two utilities are at your disposal to help with the merging:
 * `yodamerge`/`yodastack` are general-purpose scripts which work for any `YODA` file
 (but have some built-in assumptions), with usage documented below.
 [\[basic tutorial\]](merging3.md)
 * `rivet-merge` delegates merging of files back to the `Rivet` analyses which produced them.
 [\[basic tutorial\]](merging2.md)
As always, the `--help` flag will also give a lot of information about the respective script
and its limitations.

For memory-based merging of different output streams e.g. in MPI-parallelised runs
[see here](merging_mpi.md).

## Should I use `yodamerge` or `rivet-merge`?

`yodamerge` and `yodastack` are scripts shipped with `YODA`. These scripts implement a statistically-correct merging of
"live" (i.e. histogram- and profile-type) objects. However, when it comes to "inert" (i.e. estimate- or scatter-type) objects,
there are some assumptions/choices which need to be made when merging e.g. the `Scatter*D` objects:
 * should the values of each point simply be added together?
 (this assumes that each `YODA` file to be merged was generated with the same number of events)
 * should the average be taken for each point? (this assumes that each `YODA` file to be merged was generated with the same number of events)
 * perhaps the points should not be added together, but instead the list of points of all `Scatter` objects be concatenated?
 * or finally, one could even just pick the `Scatter` from the first input file and ignore the others.
The answer often depends on the details of the `finalize` method of the parent analysis.
Consider a simple efficiency (a `Scatter2D`) that is constructed from two histograms (`Histo1D` objects).
If only the resulting scatters are written out, the statistical correlations are lost and
it will be impossible to merge the files "correctly". An average might come close,
but is often not satisfactory.

This is where `rivet-merge` comes in. This script makes use of the fact that the live objects
prior to running finalize are written out to the file with the prefix `/RAW` preprended to the path.
It will merge those, then re-run finalize on the merged objects, thereby achieving the correct result,
as if all events had been processed in the same run.

As a result, `rivet-merge` should **only be used with reentrant safe routines**.
To be [reentrant-safe](merging2.md#reentrant-safety), the `finalize()` method of an analysis
should be self-consistent: everything that is required to produce the desired objects
in the output file must be booked in the initialisation phase. If you try to merge
`YODA` files from non-reentrant plugins, the script will warn you that the result
will be unpredictable.

In summary, `rivet-merge` is the more sophisticated merging tool, since it has access to the
analysis logic and can actually re-run Rivet over the merged result.
Please see the [corresponding tutorial](merging2.md) for some examples.
That said, `yodamerge` and `yodastack` are good baseline merging tools that can get you far,
and in combination with a little Python-based post-processing script, a lot is possible.
See the [corresponding tutorial](merging3.md) for some examples.

