# OpenMP-based code parallelization using ChatGPT and Github Copilot

This repository contains code for evaluating the use of ChatGPT and GitHub Copilot for OpenMP-based code parallelization. It includes nine mini-applications, each parallelized manually, using ChatGPT, and using GitHub Copilot. 

### Overview
The repository includes nine OpenMP-parallelized mini-applications that represent diverse workloads from high-performance computing (HPC). The applications test the ability of LLM-based tools to generate correct and performant parallel code.

#### The mini-applications included are:
- Feynman: Uses the Feynman-Kac algorithm to solve Poisson’s equation.
- HotSpot: A thermal simulation.
- Mandelbrot: Calculates an area of the Mandelbrot set.
- MolDyn: A molecular dynamics simulation.
- Nbody: A gravitational n-body interaction simulation.
- PiCalculation: Estimates the value of pi.
- Prime: Calculates prime numbers using a naive algorithm.
- Saxpy: Performs scalar multiplication and vector addition.
- Sgemm: Dense matrix multiplication kernel.

#### Sources
The following references provide the original benchmark applications and code repositories used in this project:
```
Che S, Boyer M, Meng J, Tarjan D, Sheaffer JW, Lee SH, et al. Rodinia: A benchmark suite for heterogeneous computing. In: 2009 IEEE international symposium on workload characterization (IISWC). Ieee; 2009. p. 44–54.
Stratton JA, Rodrigues C, Sung IJ, Obeid N, Chang LW, Anssari N, et al. Parboil: A revised benchmark suite for scientific and commercial throughput computing. Center for Reliable and High-Performance Computing. 2012;127(7.2).
Burkardt J.: John Burkardt’s home page. Accessed: 2024-04-16. https://people.sc.fsu.edu/~jburkardt/.
EPCC.: Edinburgh Parallel Computing Centre Code Repository. Accessed: 2024-04-16. https://github.com/EPCCed
Pacheco P. An introduction to parallel programming. Elsevier; 2011
AMD.: HPC training examples. Accessed: 2024-04-16. https://github.com/amd/HPCTrainingExamples/.
```

### License
This project is licensed under the GNU General Public License v3.0.
