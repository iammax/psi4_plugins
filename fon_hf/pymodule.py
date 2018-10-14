import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
from text import *
from procutil import *


def run_fon_hf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    fon_hf can be called via :py:func:`~driver.energy`.

    >>> energy('fon_hf')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Your plugin's PsiMod run sequence goes here
    PsiMod.set_global_option('BASIS', 'sto-3g')
    PsiMod.set_local_option('FON_HF', 'PRINT', 1)
    energy('scf', **kwargs)
    returnvalue = PsiMod.plugin('fon_hf.so')

    return returnvalue


# Integration with driver routines
procedures['energy']['fon_hf'] = run_fon_hf


def exampleFN():
    # Your Python code goes here
    pass
