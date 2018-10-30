from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    ext_modules=[Extension("_zeformfactor", ["_zeformfactor.cxx", "zeformfactor.cxx"])],
        include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
