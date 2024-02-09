# Migration from Rivet v3.1 to Rivet v3.2

Rivet 3.2 is a minor new release series, in which we've made some big
changes that affect the programming interface (or "API") exposed to
analysis authors.

As with the last time we made such major changes,
from version 2 to version 3, the overall effect is very positive, but
code will need to be changed. We have made the changes to all the
analyses submitted to the Rivet core library, so there's another
incentive to submit your analyses, too!

The 3.1 -> 3.2 changes particularly affect histogramming -- which
now requires the fully reworked YODA v2 -- and the configuration of
some projections. Version 3.2 also adds functionality for storing and loading
analysis-specific structure data in HDF5 format, and other new
features, but as those are not changes to existing functionality they
aren't covered in this migration guide.

- Behind the scenes a lot has changed in the histograms. However, the
  user interface for fillable histograms is still the familiar set of
  `Counter`, `Histo1D/2D`, `Profile1D/2D`, etc.  If you wish, you can now
  also create higher-dimensional histograms or histograms with discrete
  axes.

- Reference data was previously stored as point+errorbar `Scatter` types,
  which did not have proper binnings. They are now stored as `BinnedEstimate`
  types when the reference data is actually binned. The API and user-facing
  behaviour should be very similar.

- The (badly named) `BinnedHistogram` type that was basically a collection
  of histograms, themselves binned in an additional variable, has been replaced
  with the `HistoGroup`. This has the same aim and basic design, but is now
  literally implemented as a binning containing a `Histo1D` in each bin. It
  can be accessed via the predefined specialisations, `Histo1DGroupPtr` and
  `Histo2DGroupPtr`.

- Any "live" types of analysis objects left at the end of the `finalize()`
  stage will be automatically converted to "inert" estimate types, usually a
  `BinnedEstimate1D`, for direct equivalent to reference data. The `/RAW`
  pre-finalize copies will still be live types, allowing re-entrant finalization
  via `rivet-merge`. If you need to perform a custom conversion, e.g. to
  make the histogram inert without scaling by bin-width, you need to do that
  in `finalize()`, e.g. via the `barchart()` function.

- Many enums have been renamed, and unclear bool options (e.g. as arguments
  to projection constructors) have been converted to more self-documenting
  enums. Scoped enum classes are now used as standard, and are defined outside
  projection scopes, so there is now always exactly one set of `::` per enum
  value (with the implicit `Rivet::` namespace). This eliminates code smells like
  the mix of two- and three-layer enums, from both the `JetAlg` and `FastJets`
  scopes when configuring jet finders. Some renaming examples:
  - `FastJets::Algo::KT`  -> `JetAlg::KT`,
  - `JetAlg::Invisibles::NOPROMPT` -> `JetInvisibles::NONPROMPT`,
  - `JetAlg::Muons::NONE` -> `JetMuons::NONE`
  - etc.

- Tau/mu as-prompt reco bools -> `AsPromptTau/Mu` enums. In general, all
  projection-config bools are now enum classes for better
  code-readability / self-documentation.

- `DressedLeptons` projection -> `LeptonFinder`. `DressedLeptons` is now a
  `vector<DressedLepton>` alias by analogy with `Particles`, `Jets`, etc..
  The promptness bool-argument to the constructor is now a `PromptOnly::YES/NO`
  enum class.
  - Note that one of the constructors used to take a `FinalState` for photons
    followed by a `FinalState` for the bare leptons: This has now been swapped.

- `WFinder` has been removed: there are too many analysis options with
  the imperfect information induced by the invisible neutrino, so you
  should do it by hand. See the Doxygen API reference and existing
  analyses for guidance on e.g. using dressed leptons and MET via the
  `LeptonFinder`, and the new `closestMatchIndex()` metafunction used with
  the `mass()` and `mT()` functions.

- `ZFinder` has been canonically renamed to `DileptonFinder` and its
  `trackPhotons` option argument has been removed: if you do not want
  to include the photons, set the dressing dR radius to a negative
  number, conventionally `-1`. The constructor arguments have been
  changed, so the target mass now must be specified, and the
  mass-window and other cuts on the resulting dilepton vector are
  given via a single `Cut` argument. Again, see the Doxygen API
  reference and existing analyses for illustration.

- Jet smearing via `SmearedJets` now accepts multiple ordered smearing
  and efficiency functions via a parameter pack. This has required a
  re-ordering of parameters in a slightly non-optimal fashion, but was
  considered an overall preferable solution, especially given that
  this feature is a niche requirement!

- Several built-in analyses have been renamed, especially analyses for
  computing heavy-ion centrality curves, and `MC_*` decay
  analyses. These now match the `ALL_CAPS` convention used for other
  analyses, and a new conventional ordering of a `CENT` name-element
  for centrality plugins.  The few analysis base classes, e.g. for all
  `MC_*JETS` and `MC_*SPLITTINGS` analyses, have also been renamed
  into a consistent scheme.

- Any old analyses previously named with an `Sxxxxxx` publication identifier
  corresponding to the legacy SPIRES database has now been renamed to
  match the newer analyses with `Ixxxxxx` Inspire-based codes. Previous
  aliasing to allow loading via the Inspire name has now been switched so
  the non-canonical alias is the SPIRES name. These fallback aliases will
  eventually be removed.

We recommend that anyone writing a new analysis should start by running the
`rivet-mkanalysis` script rather than directly copying an old analysis code: we
update the templates generated by the script as new API features and ideas
occur, and it always represents best-practice. You can still copy-paste chunks
from your old analysis as you need them, modulo the unavoidable API changes above.
