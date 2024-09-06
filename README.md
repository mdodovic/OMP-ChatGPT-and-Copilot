# OpenMP-based code parallelization using ChatGPT and Github Copilot

This repository contains code for evaluating the use of ChatGPT and GitHub Copilot for OpenMP-based code parallelization. It includes nine mini-applications, each parallelized manually, using ChatGPT, and using GitHub Copilot. <!-- This code and data are used for the experiments presented in the paper titled: ["..."](...). -->


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
1. Che S, Boyer M, Meng J, Tarjan D, Sheaffer JW, Lee SH, et al. Rodinia: A benchmark suite for heterogeneous computing. In: 2009 IEEE international symposium on workload characterization (IISWC). Ieee; 2009. p. 44–54.
```
```
2. Stratton JA, Rodrigues C, Sung IJ, Obeid N, Chang LW, Anssari N, et al. Parboil: A revised benchmark suite for scientific and commercial throughput computing. Center for Reliable and High-Performance Computing. 2012;127(7.2).
```
```
3. Burkardt J.: John Burkardt’s home page. Accessed: 2024-04-16. https://people.sc.fsu.edu/~jburkardt/.
```
```
4. EPCC.: Edinburgh Parallel Computing Centre Code Repository. Accessed: 2024-04-16. https://github.com/EPCCed
```
```
5. Pacheco P. An introduction to parallel programming. Elsevier; 2011
```
```
6. AMD.: HPC training examples. Accessed: 2024-04-16. https://github.com/amd/HPCTrainingExamples/.
```

### Directory Structure

This repository is organized into directories that contain different versions of the mini-applications. These versions include a sequential baseline, manually parallelized code, and versions parallelized using LLMs: ChatGPT and GitHub Copilot. Below is a breakdown of the directory structure.

##### base/ 
Contains the non-parallelized, sequential version of each mini-application. These are the baseline applications used for comparison.

##### omp/
Contains manually parallelized versions of the applications using OpenMP. These implementations were optimized by human experts.

##### omp_chat_gpt/
Contains OpenMP-parallelized versions of the applications produced using ChatGPT. Only partial context (a specific code snippet) was provided to ChatGPT for these implementations.

##### omp_chat_gpt2/
Contains OpenMP-parallelized versions of the applications produced using ChatGPT with full context. The entire application code was provided to ChatGPT, allowing it to generate a more comprehensive solution.

##### copilot/
Contains OpenMP-parallelized versions of the applications produced using GitHub Copilot. The applications were parallelized based on suggestions from Copilot within Visual Studio Code.

<!-- 
### Citation
You can cite this paper as:
... 
with the following BibTeX code:
...
-->

### License
This project is licensed under the GNU General Public License v3.0.
