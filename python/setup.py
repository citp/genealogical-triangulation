from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np

# ext_modules = [Extension("cython3", ["recomb_helper.pyx",
#                                     "common_segments.pyx"])]
                         # include_dirs = [np.get_include()])

# ext_modules = [Extension("recomb_helper", ["recomb_helper.pyx"],
#                          include_dirs = [np.get_include()]),
#                Extension("common_segments", ["common_segments.pyx"],
#                          include_dirs = [np.get_include()])]

setup(
    name = 'Genetic Privacy',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    # ext_modules = ext_modules
    ext_modules = cythonize("*.pyx", include_path = [np.get_include()])
)
