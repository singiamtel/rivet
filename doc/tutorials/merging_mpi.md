
# Merging MPI-parallelised Rivet runs

When analysing events produced on multiple cores there are two common appraches
one could take to merge the Rivet output into a single output file at the end
of the run:

## File-based merging

The "brute force" approach is to initialise an `AnalysisHandler` per rank
and have each rank write out a YODA file. At the end of the runs the
various YODA files can be merged using the `rivet-merge` script
(see additional file-based merging documentation [here](merging.md)).

## Merging AnahalysisHandlers in memory

Disk space is expensive, however, and it might be more attractive to merge
the output from each individual rank in memory first, such that only a
single file needs writing out at the very end.
This can be achieved by using the `AnalysisHandler::writeData(stream, format)`
method to extract the analysis objects from the `AnalysisHandler` into a
byte stream for any given rank, which can then be sent across nodes
and gathered on a single rank. Once collected, it's straightforward to merge
the different streams using `AnalysisHandlers` together with the
`AnalysisHandler::readData(stream, format, preload)` method.
Here it is necessary to read in the complete stream by setting the `preload` flag to false
which will load every single analysis object from the stream into the `AnalysisHandler`.
The `AnalysisHandler` ojects can the be merged sequentially, followed by a
reentrant run on the merged `AnalysisHandler` in order to `finalize()` and
write out the analysis objects to an output file in the end.

In Python, this could look something like this:

```
from mpi4py import MPI
import rivet, io

def processRank(rank):
  ah = rivet.AnalysisHandler("AH%i" % rank)
  # ... analyse some events ...
  ah.collapseEventGroup()
  out = io.StringIO()
  ah.writeData(out)
  return out

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

res = processRank(mpi_rank)
res = res.getvalue().encode("utf-8")
res = mpi_comm.gather(res)

if mpi_rank == 0:
  ahmerge = rivet.AnalysisHandler("AHMERGE")
  ahmerge.readData(res[0], preload = False)
  for stream in res[1:]:
    ahtemp = rivet.AnalysisHandler("AHTEMP")
    ahtemp.readData(stream, preload = False)
    ahmerge.merge(ahtemp)
  ahmerge.finalize()
  ahmerge.writeData("mpi_merged_output.yoda.gz")
```

## (De-)Serializing AnalysisHandlers on HPCs

It's also possible to serialise the content of the `AnalysisHandler` into
a list of `double`s on each rank, which scales significantly better on HPC clusters.
Of course the deserialisation must then happen with an `AnalysisHandler` that has been
initialised with the same analyses.

In Python, this could look something like this:

```
from mpi4py import MPI
import rivet

def processRank(rank):
  ah = rivet.AnalysisHandler("AH%i" % rank)
  # ... analyse some events ...
  return ah.serializeContent()

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

res = processRank(mpi_rank)
res = mpi_comm.gather(res)

if mpi_rank == 0:
  ahmerge = rivet.AnalysisHandler("AHMERGE")
  ahtemp = rivet.AnalysisHandler("AHTEMP")
  # ... initialise with same analyses as each rank ...
  ahmerge.deserializeContent(res[0])
  xs = ahmerge.nominalCrossSection()
  xserr = ahmerge.nominalCrossSectionErr()
  for data in res[1:]:
    ahtemp.deserializeContent(data)
    xs = ahtemp.nominalCrossSection()
    xserr = ahtemp.nominalCrossSectionErr()
    ahmerge.merge(ahtemp)
  ahmerge.setCrossSection(xs, xserr)
  ahmerge.finalize()
  ahmerge.writeData("mpi_merged_output.yoda.gz")
```
