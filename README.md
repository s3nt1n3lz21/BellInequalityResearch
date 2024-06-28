# Calculation Of Tight Bell Inequalities Through Robustness

## Author: Neil Smith

**Date:** 26/07/2017  
**Contact Info:** [neilcsmith21@gmail.com](mailto:neilcsmith21@gmail.com)  
**Affiliation:** The University Of Sheffield

## Publications

### Documents
- [Thesis Part 1](./docs/thesisPart1.pdf)
- [Thesis Part 2](./docs/thesisPart2.pdf)
- [Thesis Poster](./docs/poster.pdf)

![Thesis Poster](./docs/poster.jpg)

## Description

This algorithm calculates the maximum classical value and dimension for any Bell inequality of scenarios (n, d, m), where:
- **n**: Number of parties
- **d**: Number of measurement outcomes per party
- **m**: Number of measurement settings per party

The algorithm iterates over all possible combinations of deterministic probabilities (`D_ni(di, mi)`), represented as binary variables (0 or 1). This iterative approach explores `2^(nmd)` states, making it exponential in time but polynomial in memory.

To ensure correctness and efficiency, the algorithm utilizes recursion and a class design to manage data efficiently, ensuring minimal memory usage by maintaining a single copy of the data accessible across recursive function calls.

## Usage

To use the algorithm:
1. Specify the scenario (n, d, m) and the correlator coefficient list when calling the `calcdimandclassicalbound` function.
2. The correlator coefficient list should be a row vector of length `(m+1)^n`.
3. For example, for the CHSH inequality in scenario (2, 2, 2), the list `[0 0 0 0 1 1 0 1 -1]` represents correlators:
   - `0*[00] + 0*[01] + 0*[02] + 0*[10] + 1*[11] + 1*[12] + 0*[20] + 1*[21] - 1*[22]`

### Example

An example application is provided in the `example` file, demonstrating how to apply the algorithm to specific data scenarios.

### Files

- **cbanddimcalc.m**: Algorithm implemented as a class.
- **calcdimandclassicalbound.m**: Function acting as a wrapper for the class, for user interaction.
- **main.m**: Example application demonstrating the algorithm's usage.
- **example**: Example application of the algorithm to data.
- **tests**: Test suite for validating the algorithm.

### Todos:

**Optimization:**
- Optimize array resets by clearing only necessary rows.
- Investigate using handle variables for better performance.

**Functionality:**
- Generalize for any number of measurements and outcomes per party.
- Implement a function to retrieve settings from an index using nested loops.
- Extend tensor product calculations in `main` to d-qubit states.
- Modify code to calculate either `smax` or `smax` and dimension.

**Bugs/Validation:**
- Investigate discrepancies for n = 4 compared to literature.
- Add validation to ensure correlator coefficient list isn't all zeros.
- Validate algorithm behavior for the first index.
- Improve precision by checking coefficients close to zero.