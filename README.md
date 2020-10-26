# [ARCHIVE/2008] Hartree-Fock-Roothan (F77)
  - **NAIVE** implementation of [Hartree-Fock method](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method) in [GTO base](https://en.wikipedia.org/wiki/Gaussian_orbital),
  - **ARCHIVAL** code written in **2008**,
  - pure **FORTRAN 77** syntax,
  - released as **PUBLIC DOMAIN** - use for any purpose.

# Features:
  - **NAIVE SCF** loop - O(N^4) complexity,
  - calculates **TOTAL ENERGY** and **ORBITALS** (represented by coefficients in the matrix),
  - read molecule **GEOMETRY** from file,
  - **GTO BASIS** with parameters read from file,
  - write **RESULT** to the output file.

# Limitations:
  - Only integrals for [S orbitals](https://en.wikipedia.org/wiki/Gaussian_orbital) are implemented,
  - the total number of electrons in the system **MUST BE EVEN** (see [restricted Hartree-Fock](https://en.wikipedia.org/wiki/Restricted_open-shell_Hartree%E2%80%93Fock)).

# Build:
  ```
  cd src
  make
  ```

# Usage:
  - Create input file containing:
    - molecule geometry in cartessian system,
    - basis name to use (basis must be declared within *basis.lib* file),
    - example *h2.inp* file is presented below:
  ```
GEOM
 H   0.0        0.0    0.0
 H   0.742724   0.0    0.0
END_GEOM
6-31G
```

  - Run HFR code (please not that, there is **NO** *.inp* extension in the command):
  ```
  $./hfr h2
 iter=           1    Etotal=  -1.0742384810132233
 iter=           2    Etotal=  -1.1254203943379362       delta=   5.1181913324712891E-002
 iter=           3    Etotal=  -1.1266829266059391       delta=   1.2625322680028805E-003
 iter=           4    Etotal=  -1.1267113293037336       delta=   2.8402697794582821E-005
 iter=           5    Etotal=  -1.1267119595534241       delta=   6.3024969043645740E-007
 iter=           6    Etotal=  -1.1267119735099635       delta=   1.3956539479664798E-008
```

  - Orbitals coefficients should be written to the *h2.out* file:
  ```
 Etotal=  -1.1267119735099635

 --- Funkcje falowe ---
    1.120426    0.326507    0.767506   -0.123166
   -1.347337    0.272471   -0.686424   -1.708117
   -1.120426    0.326507    0.767506    0.123166
    1.347337    0.272471   -0.686424    1.708117
```
