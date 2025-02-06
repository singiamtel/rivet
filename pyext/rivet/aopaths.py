"A collection of functions and objects for manipulating Rivet's histogram-path conventions"

def isRefPath(path):
    return path.startswith("/REF")

def isTheoryPath(path):
    return path.startswith("/THY")

def isRawPath(path):
    return path.startswith("/RAW")

def isTmpPath(path):
    "Match *any* underscore-prefixed or TMP path component"
    return "/_" in path or "/TMP/" in path


def isRefAO(ao):
    return int(ao.annotation("IsRef")) == 1 or isRefPath(ao.path())

def isRawAO(ao):
    return isRawPath(ao.path())

def isTheoryAO(ao):
    return isTheoryPath(ao.path())

def isTmpAO(ao):
    return isTmpPath(ao.path())


def stripPrefix(path):
    "Remove any ref, theory, or raw prefix dir"
    import re
    return re.sub(r"^/(REF|THY|RAW)/", "/", path)

def stripOptions(path):
    import re
    return re.sub(r":\w+=[^:/]+", "", path)

def stripWeightName(path):
    import re
    return re.sub(r"\[.*\]", "", path)

def extractWeightName(path):
    import re
    re_weight = re.compile(r".*\[(.*?)\].*")
    m = re_weight.match(path)
    if not m:
        return ""
    else:
        return m.group(1)


def extractOptionString(path):
    import re
    re_opts = re.compile(r"^.*(:\w+=[^:/]+)+")
    m = re_opts.match(path)
    if not m:
        return ""
    opts = list(m.groups())
    for i in range(len(opts)):
        opts[i] = opts[i].strip(":")
    return " [" + ",".join(opts) + "]"


class AOPath(object):
    """
    Object representation of analysis-object path structures.

    TODO: move to YODA?
    """
    import re
    re_aopath = re.compile(r"^(/[^\[\]\@\#]+)(\[[A-Za-z\d\._=\s+-]+\])?(#\d+|@[\d\.]+)?$")

    def __init__(self, path):
        self.origpath = path
        m = self.re_aopath.match(path)
        if not m:
            raise ValueError("Supplied path '%s' does not meet required structure" % path)
        self._basepath = m.group(1)
        self._varid = m.group(2).lstrip("[").rstrip("]") if m.group(2) else None
        self._binid = int(m.group(3).lstrip("#")) if m.group(3) else None
        self._isref = isRefPath(self._basepath)

    def basepath(self, keep_prefix=False):
        "Main 'Unix-like' part of the AO path, optionally retaining any prefix dir"
        p = self._basepath.rstrip("/")
        if not keep_prefix and p.startswith("/REF"):
            p = stripPrefix(p)
        return p

    def varpath(self, keep_prefix=False, defaultvarid=None):
        "The basepath, plus any bracketed variation identifier"
        p = self.basepath(keep_prefix)
        if self.varid(defaultvarid) is not None:
            p += "[%s]" % str(self.varid(defaultvarid))
        return p

    def binpath(self, keep_prefix=False, defaultbinid=None, defaultvarid=None):
        "The varpath, plus any #-prefixed bin number identifier"
        p = self.varpath(keep_prefix, defaultvarid)
        if self.binid(defaultbinid) is not None:
            p += "#%d" % self.binid(defaultbinid)
        return p

    def basepathparts(self, keep_prefix=False):
        "List of basepath components, split by forward slashes"
        return self.basepath(keep_prefix).strip("/").split("/")

    # TODO: basepathhead, basepathtail

    def dirname(self, keep_prefix=False):
        "The non-final (i.e. dir-like) part of the basepath"
        import os
        return os.path.dirname(self.basepath(keep_prefix))

    def dirnameparts(self, keep_prefix=False):
        "List of dirname components, split by forward slashes"
        return self.dirname(keep_prefix).strip("/").split("/")

    def basename(self):
        "The final (i.e. file-like) part of the basepath"
        import os
        return os.path.basename(self._basepath)

    def ananame(self):
        "The first non-prefix part of the basepath"
        parts = self.dirnameparts()
        return parts[0] if parts else None

    def varid(self, default=None):
        "The variation identifier (without brackets) if there is one, otherwise None"
        return self._varid if self._varid is not None else default

    def binid(self, default=None):
        "The bin identifier (without #) if there is one, otherwise None"
        return self._binid if self._binid is not None else default

    def isref(self):
        "Is there a /REF prefix in the original path?"
        return self._isref

    def istmp(self):
        "Do any basepath components start with an underscore, used to hide them from plotting?"
        return isTmpPath(self.basepath())

    def israw(self):
        "Do any basepath components start with /RAW, used to hide them from plotting?"
        return isRawPath(self.basepath())
