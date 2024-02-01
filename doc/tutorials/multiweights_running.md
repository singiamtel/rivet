# Run-time manipulation of weight streams

By default Rivet will run over all the available multiweights in can find in the input file,
unless they are prefixed with `IRREG_`, in which case they will always be deselected.

## Nominal weight

Rivet will run a heuristic to identify the nominal weight.
The heuristic is as follows: If it finds the empty string
among the weight names, it will treat it as the nominal weight.
If it doesn't find the empty string, it will then look for weights
labelled `WEIGHT`, `DEFAULT`, `NOMINAL` or `0`. The first match will
be treated as the nominal.

If it finds more than one, it will use the first one and print a warning
that there were potentially other candidates, thereby giving the option
to abort and re-run using a suitable `--nominal-weight` value. If it can't
find one (that's a bug in the generator HepMC interface ...), it checks if
the user has supplied a custom weight name (using the `--weight-nominal` flag)
that ought to be treated as the nominal and uses that instead.

If the user didn't supply one, Rivet falls back to the assumption that
the first weight is the nominal (technically a HepMC convention) but
will write out all weights streams in "variation weight" style
(i.e. with `[weight_name]` suffix in the path). A warning is printed.

If it can find the nominal weight stream and the user supplied a custom
one that isn't the actual nominal weight stream, Rivet will internally
continue to use the actual nominal (e.g. for the cross-section scaling)
but write out the user-supplied weight stream in "nominal weight" style
(and the actual nominal in "variation weight" style).


## Running over a subset of the available weight streams

It might be beneficial to only run on a subset of the available multiweights, especially
when the output file would be rather large and the weights aren't actually needed.
The following flags can help:

* the `--weight-match` flag can be used to select a subset of the weights,
e.g. `--weight-match=WeightName1,WeightName2`, or using a suitable regular expression,
e.g. `--weight-match=MUR.*MUF.*PDF123456`. Note that the nominal weight cannot be de-selected.
* the `--weight-unmatch` flag can be used to de-select a subset of the weights,
e.g. `--weight-unmatch=WeightName1,WeightName2`, or using a suitable regular expression,
e.g. `--weight-unmatch=MUR1_MUF1_PDF123.*`. Note that the default weight cannot be de-selected.

When both flags are passed, `--weight-match` is evaluated first, and `--weight-unmatch` is applied
to the "surviving weights".


### Skipping all variation weight streams

Passing the `--weight-skip` flag will unmatch all weights, save the nominal weight.
Only one weight stream is written out to file (in "nominal weight" style).

If the user additionally supplied `--weight-match` and `--weight-unmatch` flags,
chances are they don't really know what they're doing.
Preference is given to `--weights-skip` in this case.

If the user also supplies the `--weight-nominal` flag and asks that a variation weight
be treated as the nominal, two weight streams are written out to file: the actual default
weight stream (needed for scaling) and the custom nominal weight stream. The former will
be written out in "variation weight" style, the latter in "nominal weight" style.


