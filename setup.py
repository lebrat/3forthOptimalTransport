from distutils.core import setup, Extension
import distutils.sysconfig
import numpy

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-O2", "")

for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-DNDEBUG", "")
        
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


        
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

ot3D_module = Extension('_ot3D', libraries=['CGAL','gmp','mpfr'],
								extra_compile_args = ["-g","-std=c++11","-frounding-math","-lpthread","-lmpfr","-lgmp","-fopenmp"],
							    sources=['polyline.cxx','ot3D.cxx','LaguerreCell3D.cxx','optim.cxx','ot3D.i'],
								include_dirs = [numpy_include],swig_opts=['-c++'])

setup(name='ot3D', ext_modules=[ot3D_module], py_modules=["ot3D"])