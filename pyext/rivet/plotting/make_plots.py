import rivet, yoda
import os, glob, logging, re
from math import sqrt
from rivet.plotting import plot2yaml
from rivet.plotting.conversion_tools import type_conversion


# TODO: add more descriptive docstrings to all functions.

def _sanitise_string(s):
    s = s.replace('#','\\#')
    s = s.replace('%','\\%')
    return s


def _parse_args(args):
    """Look at the argument list and split it at colons, in order to separate
    the file names from the plotting options. Store the file names and
    file specific plotting options.

    Parameters
    ----------
    args : list[str]
        List of arguments which were previously passed to `rivet-cmphistos`.
        Format will be ['filename.yoda:key=value', ..., 'PLOT:key=value:key=value']

    Returns
    -------
    filelist : list[str]
        Raw names of the files, i.e. the first part of each string in args.
    filenames : list[str]
        Names of the files. If Name=value is passed as a plot option after a file name, this will become the filename.
        Otherwise, it will use the same value as in filelist.
    plotoptions : dict[str, dict[str, str]]
        Dictionary of plot options.
        The key will be the file name, unless Name=label is passed as a plot option after a file name,
        in which case it will be label. The dictionary value with be a dict of strings with plot options.
        One of the keys will also be PLOT (if it was passed in as an argument to args, which contains
        all plot options that will be applied to the entire figure.

    Examples
    --------
    >>> _parse_args(['mc1.yoda:Title=example title:Name=example name 1', 'mc2.yoda', 'PLOT:LogX=1'])
    (['mc1.yoda', 'mc2.yoda'],
     ['example name 1', 'mc2.yoda'],
     {
         'example name 1': {'Title': 'example title', 'Name': 'example name 1'},
         'mc2.yoda': {'Title': 'mc2'},
         'PLOT': {'LogX': '1'}
    })

    Note
    ----
    Some matplotlib line styles contain ':', which would not work with current code. TODO: change delimiter?
    """
    filelist = []
    plotoptions = {}
    for a in args:
        asplit = a.split(':')
        path = asplit[0]
        if path != "PLOT" and not re.fullmatch(r"REF(\d+)?$", path):
            filelist.append(path)
            path = path.split('/')[-1]
            while path in plotoptions:
                path = '_' + path
        plotoptions[path] = {}
        has_title = False
        has_name = ""
        for i in range(1, len(asplit)):
            ## Add 'Title' if there is no = sign before math mode
            if '=' not in asplit[i] or ('$' in asplit[i] and asplit[i].index('$') < asplit[i].index('=')):
                asplit[i] = 'Title=%s' % asplit[i]
            if asplit[i].startswith('Title='):
                has_title = True
            key, value = asplit[i].split('=', 1)
            plotoptions[path][key] = type_conversion(value)
            if asplit[i].startswith('Name=') and path != "PLOT":
                has_name = asplit[i].split('=', 1)[1]
        if has_name != "":
            plotoptions[has_name] = plotoptions.pop(path)
        else:
            has_name = path
        if path != "PLOT" and not re.fullmatch(r"REF(\d+)?$", path) and not has_title:
            plotoptions[has_name]['Title'] = _sanitise_string(os.path.basename( os.path.splitext(path)[0] ))
    return filelist, plotoptions


def _get_histos(filelist, plotoptions, path_patterns = [], path_unpatterns = [], rivetrefs = True):
    """Loop over all input files. Only use the first occurrence of any REF-histogram
    and the first occurrence in each MC file for every MC-histogram."""

    hpaths = []; anas = []
    hasVariations = False
    refhistos, mchistos = {}, {}
    for infile, inname in zip(filelist, plotoptions.keys()):
        mchistos.setdefault(inname, {})
        try:
            analysisobjects = yoda.read(infile, patterns=path_patterns, unpatterns=path_unpatterns)
        except IOError as e:
            print("File reading error:", e.strerror)
            sys.exit(1)
        for path, ao in analysisobjects.items():

            # Make a path object and ensure the path is in standard form.
            try:
                aop = rivet.AOPath(path)
            except Exception as e:
                print("Found analysis object with non-standard path structure:", path, "... skipping")
                continue

            ## We don't plot data objects with path components hidden by an underscore prefix
            if aop.istmp() or aop.israw():
                continue

            hasVariations |= bool(aop.varid())

            # Convert non-scatter objects to scatter
            ao = yoda.plotting.utils.mkPlotFriendlyScatter(ao)

            ## Add it to the ref or mc paths, if this path isn't already known
            basepath_with_anaopts = aop.basepath(keep_prefix=False)
            basepath = rivet.stripOptions(basepath_with_anaopts)
            defaultWeightName = plotoptions[inname].get('DefaultWeight', '0')
            if aop.isref() and basepath not in refhistos:
                ao.setPath(aop.varpath(keep_prefix=False, defaultvarid=defaultWeightName))
                refhistos[basepath] = ao
            else: #if basepath not in mchistos[infile]:
                mchistos[inname].setdefault(basepath_with_anaopts, {})[aop.varid(defaultWeightName)] = ao
                if basepath and basepath not in hpaths:
                    hpaths.append(basepath)
                    ana = basepath.split("/")[1]
                    if ana not in anas:
                        anas.append(ana)

    if rivetrefs:
        # Scrape Rivet ref-data files for matching ref-data AOs
        rivet_data_dirs = rivet.getAnalysisRefPaths()
        dirlist = list(set([ item for a in anas
                                  for d in rivet_data_dirs
                                  for item in glob.glob(os.path.join(d, a+'.yoda*')) ]))
        for infile in dirlist:
            analysisobjects = yoda.read(infile, patterns=path_patterns, unpatterns=path_unpatterns)
            for path, ao in analysisobjects.items():
                aop = rivet.AOPath(ao.path())
                if not aop.isref():  continue
                ao.setPath(aop.basepath(keep_prefix=False))
                new_aop = ao.path()
                if new_aop not in hpaths:  continue
                if new_aop not in refhistos:
                    refhistos[new_aop] = ao

    # propagate the variations Boolean to the top-level script
    hpaths.append(hasVariations)

    return refhistos, mchistos, hpaths


def _add_ref_hist(output, refhisto, reftitle, ratiolabel):
    refhisto.setAnnotation('IsRef', True)
    # assemble array of REF/REF2/... for every reference curve tag on
    # the command line, sort them based on tailing digits, then loop
    for ref_key in sorted(list(set([ tag for tag in output.keys()
                       if re.fullmatch(r"REF(\d+)?$", tag) ] + ['REF'])),
                       key=lambda k : int(k[3:] if len(k) > 3 else 0)):
        reflabel = reftitle + (ref_key[3:] if len(ref_key) > 3 else '')
        errpat = output.get(ref_key, {}).get('ErrorPattern', '')
        # check if a band was explicitly requested
        needsBand = len([ val for key, val in output.get(ref_key, {}).items()
                                           if 'ErrorBand' in key and val ])
        # add new reference curve
        output['histograms'][reflabel] = {
            'nominal' : yoda.plotting.utils.mkPlotFriendlyScatter(refhisto, errorPattern=errpat)
        }
        output['histograms'][reflabel]['IsRef'] = True
        output['histograms'][reflabel]['LineColor'] = 'black'
        # set label for reference data in legend
        output['histograms'][reflabel]['Title'] = reftitle
        # set additional reference-data options
        output['histograms'][reflabel].update(output.get(ref_key, {}))
        if needsBand:
            bandao = yoda.plotting.utils.mkPlotFriendlyScatter(refhisto, errorPattern=errpat)
            output['histograms'][reflabel]['BandUncertainty'] = bandao

    # decide if ratio panel is shown or not
    output['plot features']['RatioPlot'] = ratiolabel != None
    # set label on y-axis of the ratio panel
    output['plot features']['RatioPlotYLabel'] = ratiolabel


def _make_output(plot_id, plotdirs, config_files, mchistos, refhistos, plotoptions,
                 style, rc_params, mc_errs, nRatioTicks, showWeights, removeOptions, deviation,
                 canvasText, refLabel=None, ratioPlotLabel = None, showRatio = None, verbose = False,):

    """Create output dictionary for the plot_id.

    Parameters
    ----------
    plot_id : str
        ID, usually of the format AnalysisID/HistogramID.
    plotdirs : list[str]
        All directories to look for .plot files at.
    config_files : list[str]
        Additional plot settings that will be applied to all figures.
    mchistos : dict
        Dictionary of the Monte Carlo YODA histograms.
        The structure is {filename: {plot_id: {"0": yoda_histogram1, "1": yoda_histogram2, ...}}}
        Usually only "0" exists as the innermost key.
    refhistos : dict
        Dictionary of the reference analysis data YODA histograms.
    plotoptions : dict[str, dict[str, str]]
        Dict containing all plot options for all histograms and all plots.
    mc_errs : bool
        See rivet_mkdat
    style : str
        A predefined name of a style.
    removeOptions : bool
        If true, prevents appending the options string to the legend label
    deviation : bool
        If true, express compatability between curve and ref. data in terms
        of standard deviations in ratio panel.
    rc_params : dict[str, str]
        Dict of rcParams that will be added to the rcParams section of the output .dat file.

    Returns
    -------
    outputdict : dict
        Correctly formatted dictionary that can be passed to `yaml.dump` to write to an output file.
    """
    outputdict = {}
    plot_configs = plot2yaml.get_plot_configs(plot_id, plotdirs=plotdirs, config_files=config_files)
    outputdict['plot features'] = plot_configs
    if deviation:
        outputdict['plot features']['RatioPlotMode'] = 'deviation'
    rpmode = plot_configs.get('RatioPlotMode', 'mcdata')

    # only write extra info to the .dat file if specified by user
    if nRatioTicks !=1: outputdict['plot features'].update({"nRatioTicks": nRatioTicks})
    if canvasText != None: outputdict['plot features'].update({"canvasText" : canvasText})
    outputdict['plot features'].update(plotoptions.get('PLOT', {}))
    outputdict['rcParams'] = rc_params
    outputdict['style'] = style
    outputdict['stylepath'] = '../'
    outputdict['histograms'] = {}
    outputdict['pat_warn'] = {}
    for key, val in plotoptions.items():
        if re.fullmatch(r"REF(\d+)?$", key):
            outputdict[key] = val

    componentNames = ['BandComponentPDF', 'BandComponentEnv']

    # Check if there's reference data
    if plot_id in refhistos and rpmode != 'datamc':
        # Adding it first ensures it goes in the denominator
        rplabel = None
        hasRatio = showRatio if showRatio != None else refhistos[plot_id].annotation('RatioPlot', True)
        if hasRatio:
            rplabel = ratioPlotLabel if ratioPlotLabel != None else \
                      refhistos[plot_id].annotation('RatioPlotYLabel', 'MC/Data')
        reftitle = refLabel if refLabel != None else \
                   refhistos[plot_id].annotation('Title', 'Data')
        reflabel = reftitle if reftitle != None else 'Data'
        _add_ref_hist(outputdict, refhistos[plot_id], reflabel, rplabel)

    # Now add MC curves
    lhapdfCheck = True
    for filename, mchistos_in_file in mchistos.items():
        for plot_id_with_anaopt in sorted(mchistos_in_file):
            if rivet.stripOptions(plot_id_with_anaopt) != plot_id:
                continue
            histogroup = mchistos_in_file[plot_id_with_anaopt]

            label = rivet.extractOptionString(plot_id_with_anaopt)
            outputdict['histograms'][filename+label] = {}

            thisFilePlotOptions = dict(plotoptions.get(filename, {}))
            # add options string to legend entry
            newtitle = thisFilePlotOptions.get('Title', '')
            if not removeOptions:
                newtitle += label
            thisFilePlotOptions['Title'] = newtitle
            outputdict['histograms'][filename+label].update(thisFilePlotOptions)

            makePDFBand  = thisFilePlotOptions.get('BandComponentPDF', '')
            makeEnvelope = thisFilePlotOptions.get('BandComponentEnv', '')

            # check if lhapdf is available
            if makePDFBand and lhapdfCheck:
                try:
                    import lhapdf
                    lhapdf.setVerbosity(0)
                    lhapdfCheck = False
                except ImportError as e:
                    print("LHAPDF not available! Need this to construct PDF band:",
                            f" failing `import {e.name}`")
                    exit(1)

            nominalVariationKey = thisFilePlotOptions.get('DefaultWeight', '0')
            if nominalVariationKey == None:
                raise NameError("Could not find nominal variation weight!")

            nomVals = None
            pdf_matches = { }; env_matches = { }
            PDFvars = [ [] for _ in makePDFBand.split() ]
            PDFsets = [ None for _ in makePDFBand.split() ]
            Enverrors = [ [] for _ in makeEnvelope.split() ]
            for histogramkey, histogram in histogroup.items():

                if histogramkey.startswith('EXTRA'):  continue  # cf. weightname convention

                isNominal = (nominalVariationKey == histogramkey)
                # Maybe add this mc_errs option to the plotoptions dict and only
                # pass the plotoptions dict to the function?
                outputdict['histograms'][filename+label]['ErrorBars'] = mc_errs

                thisObj = yoda.plotting.utils.mkPlotFriendlyScatter(histogram)

                # no support for 3D scatters and bands
                if thisObj.type() == "Scatter3D":
                    makeEnvelope = makePDFBand = ''

                # central values of current object
                central_values = thisObj.vals(thisObj.dim()-1)

                if isNominal:
                    nominalScatter = thisObj
                    outputdict['histograms'][filename+label]['nominal'] = thisObj

                for i, prescription in enumerate(makePDFBand.split()):
                    if prescription not in pdf_matches:
                        pdf_matches[prescription] = [ ]
                    if any([ re.search(pat, histogramkey) for pat in prescription.split(',') ]):
                        pdf_matches[prescription].append(histogramkey)

                        # store values from PDF variation
                        PDFvars[i].append(central_values)

                        # initialize pdf set object from lhapdf
                        if PDFsets[i] is None:
                            lhapdfID = int(re.search('PDF[0-9]*', histogramkey).group(0)[3:])
                            PDFsets[i] = lhapdf.mkPDF(lhapdfID).set()

                for i, prescription in enumerate(makeEnvelope.split()):
                    if prescription not in env_matches:
                        env_matches[prescription] = [ ]
                    if isNominal or any([ re.search(pat, histogramkey) for pat in prescription.split(',') ]):
                        if not isNominal:
                            env_matches[prescription].append(histogramkey)
                        if not Enverrors[i]:
                            Enverrors[i] = [ list(central_values), list(central_values) ]
                        else:
                            Enverrors[i][0] = list(map(min, zip(Enverrors[i][0], central_values)))
                            Enverrors[i][1] = list(map(max, zip(Enverrors[i][1], central_values)))

                # Don't plot multiweights if already plotting a band
                if showWeights and not isNominal and not makeEnvelope and not makePDFBand:
                    # Check if the user supplied regex-based weightname filtering
                    var_filter = thisFilePlotOptions.get('Variations', '')
                    if var_filter:
                        if not any([ re.search(pat, histogramkey) for pat in var_filter.split(',') ]):
                            continue
                    outputdict['histograms'][filename+label]['multiweight'+histogramkey] = thisObj

            for pat in pdf_matches:
                if not len(pdf_matches[pat]):
                    if verbose:
                        print (f"WARNING: PDF band prescription '{pat}' did not match any variation weights!")
                    elif not pat in outputdict['pat_warn'].setdefault('PDF band', {}).setdefault(filename, []):
                        outputdict['pat_warn']['PDF band'][filename].append(pat)
                elif verbose:
                    print ("PDF prescription \"%s\" matches:" % pat)
                    print (pdf_matches[pat])
            for pat in env_matches:
                if not len(env_matches[pat]):
                    if verbose:
                        print (f"WARNING: Envelope prescription '{pat}' did not match any variation weights!")
                    elif not pat in outputdict['pat_warn'].setdefault('Envelope', {}).setdefault(filename, []):
                        outputdict['pat_warn']['Envelope'][filename].append(pat)
                elif verbose:
                    print ("Envelope prescription \"%s\" matches:" % pat)
                    print (env_matches[pat])
            del pdf_matches, env_matches

            PDFerrors = [ ]
            for pdf_set, pdf_vars in zip(PDFsets, PDFvars):
                # if number of PDF variations if off by 1,
                # probably needs the nominal
                if len(pdf_vars) == int(pdf_set.size) - 1:
                    pdf_vars.append(nomVals)
                elif len(pdf_vars) != int(pdf_set.size):
                    raise ValueError("Number of matched PDF variations is %i, expected %s!" % (len(pdf_vars), pdf_set.size))
                pdf_vars = list(map(list,zip(*pdf_vars))) # transpose
                try:
                    # calculate uncertainties from all PDFs multiweight
                    # histos that matched regex from the user
                    uncertainties = [ pdf_set.uncertainty(binVars) for binVars in pdf_vars ]
                    PDFerrors.append([ (unc.errminus, unc.errplus) for unc in uncertainties ])

                except RuntimeError:
                    print("Error in constructing the PDFset. Skipping.")

            # let user ask for a band, if no BandComponentEnv/PDF provided
            # this will just be a band with stat. errors
            if makePDFBand or makeEnvelope:

                BandScatter = nominalScatter.clone()
                # iterate over bins
                for ibin, y in enumerate(nominalScatter.vals(nominalScatter.dim()-1)):
                    totErrDn, totErrUp = BandScatter.point(ibin).errs(BandScatter.dim()-1)
                    totErrDn = totErrDn*totErrDn
                    totErrUp = totErrUp*totErrUp

                    # add PDF uncertainty in quadrature
                    for errs in PDFerrors:
                        totErrDn += errs[ibin][0]*errs[ibin][0]
                        totErrUp += errs[ibin][1]*errs[ibin][1]

                    # add Envelope uncertainty in quadrature
                    for errDn, errUp in Enverrors:
                        absEnvDn = y - errDn[ibin]
                        absEnvUp = errUp[ibin] - y
                        totErrDn += absEnvDn*absEnvDn
                        totErrUp += absEnvUp*absEnvUp

                    # Scatter object with total Band uncertainty
                    BandScatter.point(ibin).setErrs(BandScatter.dim()-1, sqrt(totErrDn), sqrt(totErrUp))
                outputdict['histograms'][filename+label]['BandUncertainty'] = BandScatter

    # Check if there's reference data
    if plot_id in refhistos and rpmode == 'datamc':
        # Adding it last ensures it goes in the numerator
        rplabel = None
        hasRatio = showRatio if showRatio != None else refhistos[plot_id].annotation('RatioPlot', True)
        if hasRatio:
            rplabel = ratioPlotLabel if ratioPlotLabel != None else \
                      refhistos[plot_id].annotation('RatioPlotYLabel', 'Data/MC')
        reftitle = refLabel if refLabel != None else \
                   refhistos[plot_id].annotation('Title', 'Data')
        reflabel = reftitle if reftitle != None else 'Data'
        _add_ref_hist(outputdict, refhistos[plot_id], reflabel, rplabel)

    # Remove all sections of the output_dict that do not contain any information.
    # A list of keys is first created. Otherwise, it will raise an error since the size of the dict changes.
    dict_keys = list(outputdict.keys())
    for key in dict_keys:
        if not outputdict[key] or re.fullmatch(r"REF(\d+)?$", key):
            del outputdict[key]
    return outputdict


def assemble_plotting_data(args, path_pwd=True, rivetrefs=True,
                           path_patterns=[], path_unpatterns=[],
                           plotinfodirs=[], style='default', config_files=[],
                           hier_output=False, outdir='.', mc_errs=True,
                           rivetplotpaths=True, analysispaths=[], verbose=False,
                           nRatioTicks=1, showWeights=False,
                           removeOptions = False, deviation=False,
                           canvasText=None, refLabel=None, ratioPlotLabel=None,
                           showRatio=None):
    """Create a dictionary of the plotting data that can be turned
    into self-consistent Python executables.

    Parameters
    ----------
    args : Iterable[str]
        Non-keyword arguments that were previously passed to rivet-cmphistos.
        E.g., ['mc1.yoda', 'mc2.yoda:Title=example title', 'PLOT:LogX=1']
    path_pwd : bool
        Search for plot files and reference data files in current directory.
    rivetrefs : bool
        If False, don't use Rivet reference data files
    path_patterns : Iterable[str]
        Only write out histograms whose $path/$name string matches these regexes.
        The argument may also be a text file.
    path_unpatterns : Iterable[str]
        Exclude histograms whose $path/$name string matches these regexes
    plotinfodirs : list[str]
        Directory which may contain plot header information (in addition to standard Rivet search paths).
    style : str
        Set the style of all plots and additional rcParams.
        Format is style:key=value:key2=value2...
        The first part of the string must be a name of a builtin style (e.g. 'default').
        The other keys and values must be valid rcParams.
        However, the validity is not checked by this function.
    config_files : list[str]
        Additional plot config file(s).
        Settings will be included in the output configuration.
        ~/.make-plots will automatically be added.
    hier_output : bool
        Write output .dat files into a directory hierarchy which matches the analysis paths.
    outdir : str
        Write dat files into this directory.
    mc_errs : bool
        If True, add the errors of the Monte-Carlo histograms.
    rivetplotpaths : bool
        Search for .plot files in the standard Rivet plot paths.
    verbose : bool
        If True, write more information to stdout.
    nRatioTicks: int
        Number of minor ticks between major ticks, can be specified in rivet-mkhtml
    deviation: bool
        Scale ratio-plot to error of the reference histogram (1 standard deviation)
    refLabel : str
        Legend name of the reference data in the plots.
    ratioPlotLabel : str
        Label on the y-axis of the ratio panel.

    Returns
    -------
    dict[str, dict]
        A dict containing all dicts that are usually written to the dat file. The key is the analysis ID.

    Raises
    ------
    IOError
        If the program does not have read access to .plot or .yoda files, or if it cannot write the output .dat files.

    Notes
    -----
    TODO The keys in the returned dict always includes a / rather than being the actual output file name.
    The keys will therefore differ from the actual output file names when hier_output == False.
    To get the actual file names, / should be replaced by _ when hier_output == False.
    """

    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    # TODO: more elegant solution for getting rc_params by refactoring _parse_args.
    #  Then the 4 lines below can be replaced by 1 line
    stylename, rc_params_dict = _parse_args([style])
    stylename = stylename[0]    # Convert list to str
    rc_params_dict = rc_params_dict[stylename]  # Convert dict of dicts to dict
    del rc_params_dict['Title']

    ## Add pwd to search paths
    if path_pwd:
        rivet.addAnalysisLibPath(os.path.abspath("."))
        rivet.addAnalysisDataPath(os.path.abspath("."))
    for path in analysispaths:
        rivet.addAnalysisLibPath(os.path.abspath(path))
        rivet.addAnalysisDataPath(os.path.abspath(path))

    # Split the input file names and the associated plotting options given on the command line into two separate lists
    filelist, plotoptions = _parse_args(args)

    ## Check that the files exist
    for f in filelist:
        if not os.access(f, os.R_OK):
            raise IOError("Error: cannot read from %s" % f)

    plotdirs = plotinfodirs
    plotdirs += [os.path.abspath(os.path.dirname(f)) for f in filelist]
    plotdirs += (rivet.getAnalysisPlotPaths() if rivetplotpaths else [])

    # Create a list of all histograms to be plotted, and identify if they are 2D histos (which need special plotting)
    refhistos, mchistos, hpaths = _get_histos(filelist, plotoptions, path_patterns, path_unpatterns, rivetrefs)
    hasVariations = hpaths.pop()

    # Write each file
    plot_info_dicts, pattern_warnings = {}, {}
    for plot_id in hpaths:
        outputdict = _make_output(
            plot_id, plotdirs, config_files,
            mchistos, refhistos,
            plotoptions, stylename, rc_params_dict, mc_errs,
            nRatioTicks, showWeights, removeOptions, deviation,
            canvasText, refLabel, ratioPlotLabel, showRatio, verbose
        )
        for presc, warnings in outputdict.setdefault('pat_warn', {}).items():
            for filename, pats in warnings.items():
                for pat in pats:
                    if not pat in pattern_warnings.setdefault(presc, {}).setdefault(filename, []):
                        pattern_warnings[presc][filename].append(pat)
        del outputdict['pat_warn']
        plot_info_dicts[plot_id] = outputdict
    if not verbose:
        for presc, warnings in pattern_warnings.items():
            for filename, pats in warnings.items():
                for pat in pats:
                    print (f"WARNING: {presc} prescription '{pat}' did not match "
                           + f"any variation weights in file '{filename}'!")

    return plot_info_dicts, hasVariations

