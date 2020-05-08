import shutil 
import os

def find_executable(names, default = None):
    """
    Return the path to the first executable found in PATH from the given list
    of names.

    Raises a (hopefully helpful) error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of %s in PATH=%s" % (names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception

    return exe

find_executable("raxmlHPC-PTHREADS-AVX2")
