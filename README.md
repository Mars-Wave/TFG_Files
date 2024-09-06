# K-Anonymity Problem Solvers

This repository contains the files and scripts for my end-of-degree project, focusing on solving the **K-Anonymity** problem using multiple algorithmic approaches. The repository includes implementations of a **Genetic Algorithm**, an **Exhaustive Branch-and-Bound Algorithm**, and **Minizinc models**. Additionally, Python scripts are provided to generate plots for visualizing the results of these algorithms.

## K-Anonymity brief description
The K-Anonymity problem is a fundamental issue in data privacy, aiming to protect sensitive information by generalizing records so that they are indistinguishable from at least \(k-1\) others. This repository explores various approaches to solve the problem in the context of medical databases, at the point of the screening process. The following algorithms were designed and can be found in this repository:

## Solutions designed
- **Steady-State Genetic Algorithm**: A heuristic-based method that evolves solutions over time using techniques such as heuristic-guided mutation, crossover and gene repair.
- **Exhaustive Branch-and-Bound Algorithm**: A precise, deterministic method that explores all possible solutions.
- **Minizinc Files**: Constraint models implemented in the Minizinc modeling language to solve the problem using solvers like Gecode.
