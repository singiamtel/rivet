
# Merging MPI-parallelised Rivet runs

When analysing events produced on multiple cores there are two common appraches
one could take to merge the Rivet output into a single output file at the end
of the run:

## File-based merging

The "brute force" approach is to initialise an `AnalysisHandler` per rank
and have each rank write out a YODA file. At the end of the runs the
various YODA files can be merged using the `rivet-merge` script
(see additional file-based merging documentation [here](merging.md)).

## Merging AnalysisHandlers in memory

Disk space is expensive, however, and it might be more attractive to merge
the output from each individual rank in memory first, such that only a
single file needs writing out at the very end.

This can be achieved by using YODA's (de-)serialisation methods, which
allows representing the numerical content of the `AnalysisHandler`
as a long list of floating-point values. The streams of values can then be collapsed
using an MPI-reduce operation and deserialised into the `AnalysisHandler`
of the root/master rank. Of course the deserialisation must then happen with an
`AnalysisHandler` that has been

In Python, this could look something like this:

```
from mpi4py import MPI
import rivet, io

def processRank(rank):
  ah = rivet.AnalysisHandler("AH%i" % rank)
  # ... analyse some events ...
  ah.collapseEventGroup()
  return ah.serializeContent(True)

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

res = processRank(mpi_rank)
res = mpi_comm.reduce(res)

if mpi_rank == 0:
  ah.deserializeContent(red[0], mpi_size)
  ah.finalize()
  ah.writeData("mpi_merged_output.yoda.gz")
```

