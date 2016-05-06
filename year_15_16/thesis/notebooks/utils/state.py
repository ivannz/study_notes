import time
import pickle
import gzip

def _save(obj, path, gz=None):
    """Pickles a pythonic `obj` into a file given by `path`.

    Parameters
    ----------
    obj: any python object
        An object to pickle.
    path: string
        A file in which to pickle the object.
    gz: integer, or None
        If None, then does not apply compression while pickling. Otherwise must 
        be an integer 0-9 which determines the level of GZip compression: the lower
        the level the less thorough but the more faster the compression is.
        Value `0` produces a GZip archive with no compression whatsoever, whereas
        the value of `9` produces the most compressed archive.

    Returns
    -------
    filename: string
        The name of the resulting archive.
    """
    if not(gz is None or (isinstance(gz, int) and 0 <= gz <= 9)):
        raise TypeError("""`gz` parameter must be either `None` """
                        """or an integer 0-9.""")

    open_ = open if gz is None else lambda f, m: gzip.open(f, m, gz)
    filename_ = "%s-%s.%s"%(path, time.strftime("%Y%m%d_%H%M%S"),
                            "pic" if gz is None else "gz")

    with open_(filename_, "wb+") as f:
        pickle.dump(obj, f)
    return filename_

def _load(filename):
    """Recover an object from the file identified by `filename`.
    Automatically handles plain and GZipped data blobs.

    Parameters
    ----------
    filename: string
        A `file` in which an object is pickled.

    Returns
    -------
    object: a python object
        The recovered pythonic object.
    """
    open_ = open if not filename.endswith(".gz") else gzip.open

    with open_(filename, "rb") as f:
        obj = pickle.load(f)
    return obj
