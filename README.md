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
When the task is complete, energy consumption can be obtained with `eacct`: `eacct -j <job_id>`.

## How to cite
Our papers: 
- [Enabling mixed-precision with the help of tools: A Nekbone case study](https://link.springer.com/chapter/10.1007/978-3-031-85697-6_3)
- [Enabling mixed-precision in spectral element codes](https://www.sciencedirect.com/science/article/pii/S0167739X25002857)

Cite:

- Chen, Y., Castro, P. D. O., Bientinesi, P., & Iakymchuk, R. (2024, September). Enabling mixed-precision with the help of tools: A Nekbone case study. In International Conference on Parallel Processing and Applied Mathematics (pp. 34-50). Cham: Springer Nature Switzerland.
- Chen, Y., de Oliveira Castro, P., Bientinesi, P., Jansson, N., & Iakymchuk, R. (2025). Enabling mixed-precision in spectral element codes. Future Generation Computer Systems, 107990.

BibTeX:
```plain
@inproceedings{chen2024enabling,
title = {Enabling mixed-precision with the help of tools: A Nekbone case study},
author = {Yanxiang Chen and Pablo {de Oliveira Castro} and Paolo Bientinesi and Roman Iakymchuk}
booktitle = {International Conference on Parallel Processing and Applied Mathematics},
pages = {34--50},
year = {2024},
organization = {Springer},
doi = {https://doi.org/10.1007/978-3-031-85697-6_3},
url = {https://link.springer.com/chapter/10.1007/978-3-031-85697-6_3}
}
```

```plain
@article{CHEN2026107990,
title = {Enabling mixed-precision in spectral element codes},
author = {Yanxiang Chen and Pablo {de Oliveira Castro} and Paolo Bientinesi and Niclas Jansson and Roman Iakymchuk},
journal = {Future Generation Computer Systems},
volume = {174},
pages = {107990},
year = {2026},
issn = {0167-739X},
doi = {https://doi.org/10.1016/j.future.2025.107990},
url = {https://www.sciencedirect.com/science/article/pii/S0167739X25002857}
```
