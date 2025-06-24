# Nekbone-mp

## Supported compilers

To compile mixed-precision Nekbone, both Fortran compiler and C compiler is required. Tested compilers include:

- Classic Flang: flang & clang (https://github.com/flang-compiler/flang)
- Intel compiler: ifort & icc (our environment)
- Intel oneAPI compiler: ifx & icx
- Cray compiler: ftn & cc
- GNU compiler: gfortran & gcc (old version like gfortran-5)

## Test cases

- Test case w/o preconditioner: [test/no_precond](https://github.com/yanxchen/enabling-mixed-precision/tree/main/Nekbone-mp/test/no_precond)
- Test case w preconditioner: [test/precond](https://github.com/yanxchen/enabling-mixed-precision/tree/main/Nekbone-mp/test/precond)

## Steps:

1. Specify compilers in makenek file: `F77` and `CC`, like ifort/icc.
2. (Optional) For parallel building: comment `IFMPI="false"` and specify MPI compiler.
3. (Optional) For changing optimization flags: modify `OPT_FLAGS_STD` and `OPT_FLAGS_MAG`.
4. Run the script: `./makenek`.
5. Define the parameters for simulation in `data.rea`.
6. Run the program `./nekbone`.
