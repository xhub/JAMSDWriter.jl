[![Build Status](https://travis-ci.com/xhub/JAMSDWriter.jl.svg?branch=master)](https://travis-ci.com/xhub/JAMSDWriter.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/w4wrpfi0x2k8a7ih?svg=true)](https://ci.appveyor.com/project/xhub/jamsdwriter-jl)
[![Coverage Status](https://coveralls.io/repos/github/xhub/JAMSDWriter.jl/badge.svg?branch=master)](https://coveralls.io/github/xhub/JAMSDWriter.jl?branch=master)
[![codecov](https://codecov.io/gh/xhub/JAMSDWriter.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xhub/JAMSDWriter.jl)


# JAMSDWriter

This package enables to export a MathProgBase model to the JAMSD reformulation solver.

This package is heavily based on AmplNLWriter.jl

## Installing

This package requires the binary library jamsd to be available on the system and a valid GAMS install, with the gams executable in the system path.
For the latter, the user can refer to the GAMS installation instructions. For the former, this can be achieved in two ways:
- put the library file in the same directory as the libraries bundled with julia (usually lib/julia in the base directory of the julia install).
- modify either LD_LIBRARY_PATH (linux) or DYLD_LIBRARY_PATH (Mac OS) or PATH (Windows) to include the directory where the library is located

Note that the JAMSD solver that comes with GAMS cannot be used as for that purpose.

An optional dependency of JAMSD is the library vrepr. The latter is used to compute a vertex representation of a polyhedron, and is required
in order to use the epigraph reformulation of an OVF function. This library is dynamically loaded is this kind of computation is requested.
Follow the same install procedure as the JAMSD library.
