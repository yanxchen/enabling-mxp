# Enabling mixed-precision in spectral element codes

This repository contains source codes and softwares used in our paper: 

_Enabling mixed-precision in spectral element codes_.

## Methodology / workflow for enabling mixed-precision
<img src="https://github.com/yanxchen/enabling-mixed-precision/blob/main/methodology.png" width="550">

Arithmetic tool for accuracy check:
- Verificarlo: https://github.com/verificarlo/verificarlo

## Case studies
Please check the subfolders for detailed build steps.

### Mixed-precision Nekbone
> Nekbone solves a standard Poisson equation using a conjugate gradient iteration with a simple or spectral element multigrid preconditioner on a block or linear geometry. It exposes the principal computational kernel to reveal the essential elements of the algorithmic-architectural coupling that is pertinent to Nek5000.
> Original repository: https://github.com/Nek5000/Nekbone

[Nekbone-mp](https://github.com/yanxchen/enabling-mixed-precision/tree/main/Nekbone-mp) is an implementation of mixed-precision Nekbone, for the case w and w/o preconditioner.

### Mixed-precision Neko
> Neko is a portable framework for high-order spectral element flow simulations. Written in modern Fortran, Neko adopts an object-oriented approach, allowing multi-tier abstractions of the solver stack and facilitating various hardware backends ranging from general-purpose processors, CUDA and HIP enabled accelerators to SX-Aurora vector processors. Neko has its roots in the spectral element code Nek5000 from UChicago/ANL, from where many of the namings, code structure and numerical methods are adopted.
> Original repository: https://github.com/ExtremeFLOW/neko

[Neko-mp](https://github.com/yanxchen/enabling-mixed-precision/tree/main/Neko-mp) is an implementation of mixed-precision Neko, with the Poission equation problem as example.

## Results

### Time-to-solution measurements
We used ```time``` and timing functions inside the codes to measure time-to-solution.

### Energy-to-solution measurements
On [LUMI](https://www.lumi-supercomputer.eu/lumi_supercomputer/) supercomputer, we obtained energy data with the help of SLURM tools:
```shell
sacct --format=ConsumedEnergy -j <jobId>
```

On [MareNostrum5](https://www.bsc.es/supportkc/docs/MareNostrum5/intro/) supercompuer, we used [EAR](https://www.bsc.es/research-and-development/software-and-apps/software-list/ear-energy-management-framework-hpc) software for energy measurements.
A quick and easy approch is to use SLURM script:
```bash
#SBATCH -A <account_name>
#SBATCH -q <queue_name>
#SBATCH --job-name=Nekbone-MP_ENG
#SBATCH --exclusive                     # exclusive access
#SBATCH --nodes=1
#SBATCH --tasks-per-node=80
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=%j_Nekbone-MP_ENG.out
#SBATCH --ear=on                        # turn on EAR plugin
#SBATCH --ear-policy=monitoring         # only monitor

srun ./nekbone
```
