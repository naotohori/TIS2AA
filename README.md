[![DOI](https://zenodo.org/badge/85434321.svg)](https://zenodo.org/badge/latestdoi/85434321)

# TIS2AA
----

#### Requirements:
 * Numpy
 * Fortran compiler
 * For minimization, AmberTools (https://ambermd.org/AmberTools.php) is required.

## Getting started

### 1. Clone the repository including submodules.

```sh
$ git clone --recurse-submodule https://github.com/naotohori/TIS2AA 
```

If you already have this repository without submodules, use the following command inside the directory.

```sh
$ git submodule update --init
```

### 2. Library preparation.

```sh
$ cd TIS2AA/bin/fQCP
$ make f2py
$ cd ../..
```

A shared object `CalcROT.cpython-XX-XXXX.so` will be created. 

### 3. Download fragment libraries

```sh
$ curl -L https://www.dropbox.com/s/uuiif7di6vg8wex/RNA09_FRAG.tar.gz | tar zx
$ curl -L https://www.dropbox.com/s/5mxp18qsn568p23/DNA_FRAG.tar.gz | tar zx
```

## Usage

````bash
 $ tis2aa.py  trna.cg.pdb  trna.log trna.aa.pdb
 $ min.py  trna.aa.pdb  trna.aa.min.pdb
````

----
#### References
 1. Murray, L. J. W., Arendall, W. B., Richardson, D. C., & Richardson, J. S. (2003). RNA backbone is rotameric. Proceedings of the National Academy of Sciences, 100(24), 13904–13909. http://doi.org/10.1073/pnas.1835769100
 2. Humphris-Narayanan, E., & Pyle, A. M. (2012). Discrete RNA Libraries from Pseudo-Torsional Space. Journal of Molecular Biology, 421(1), 6–26. http://doi.org/10.1016/j.jmb.2012.03.002
