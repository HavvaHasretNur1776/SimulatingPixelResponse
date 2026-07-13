# Pixel Detector Response Simulation

A modular C++ framework for detector response simulation, digitisation, and fast track reconstruction developed within the CERN LHCb software framework.

Although implemented using the Gaudi framework, the software architecture and numerical components are designed around transferable simulation concepts including detector response modelling, numerical optimisation, efficient memory management, and modular reconstruction algorithms.

---

## Overview

This project models the response of a silicon pixel detector, transforming Monte Carlo particle hits into realistic detector measurements through a complete digitisation pipeline.

The repository contains both the detector simulation chain and a collection of reusable numerical components developed to efficiently reconstruct particle trajectories.

Key features include:

- Detector response simulation
- Pixel digitisation
- Charge deposition modelling
- Fake hit generation
- Fast polynomial-based track fitting
- Performance monitoring and validation
- Reusable numerical fitting utilities

---

## Simulation Pipeline

```
Monte Carlo Hits
        │
        ▼
Charge Deposition
        │
        ▼
Pixel Digitisation
        │
        ▼
Noise / Fake Hit Generation
        │
        ▼
Fast Track Reconstruction
        │
        ▼
Validation & Monitoring
```

---

## Repository Structure

```
include/
```

Reusable numerical components

- `PolyFitter.h` – weighted polynomial least-squares fitting
- `PolyEvaluator.h` – polynomial evaluation
- `FitResult.h` – stores fit parameters and covariance information
- `PackedSymmetricView.h` – memory-efficient symmetric matrix storage
- `myCholeskyDecomp.h` – Cholesky decomposition wrapper
- `span.h` – lightweight span implementation

```
src/
```

Detector simulation pipeline

- Monte Carlo hit processing
- Charge deposition
- Pixel digitisation
- Fake hit generation
- Track fitting
- Monitoring algorithms

---

## Numerical Components

The repository contains several reusable numerical utilities developed to minimise computational overhead during reconstruction.

### PackedSymmetricView

A lightweight matrix view storing only the lower triangular part of symmetric matrices.

Features

- approximately 50% lower memory usage compared to dense storage
- transparent matrix indexing
- zero-copy view
- cache-friendly contiguous storage

---

### PolyFitter

A generic weighted polynomial least-squares fitting framework.

Features

- template-based implementation
- incremental accumulation of normal equations
- single contiguous memory allocation
- weighted measurements
- Cholesky-based solver
- reusable for arbitrary polynomial fitting problems

---

### FitResult

Stores

- fitted coefficients
- covariance information
- polynomial evaluator
- fit metadata

---

## Performance Considerations

The fitting framework was designed with computational efficiency in mind.

Optimisations include

- packed symmetric matrix storage
- contiguous memory allocation
- reusable work buffers
- template-based generic implementation
- incremental construction of normal equations
- Cholesky decomposition for efficient linear system solving

These design choices reduce memory allocations and improve cache locality during repeated fits.

---

## Preliminary

- Modern C++
- Gaudi Framework
- ROOT
- Numerical Linear Algebra
- Cholesky Decomposition


---

## Skills Demonstrated

- Object-Oriented Design
- Numerical Optimisation
- Detector Response Modelling
- Fast Reconstruction Algorithms
- Memory-Efficient Data Structures
- Template Programming

---

## Future Work

Possible future extensions include

- GPU-based reconstruction
- Exploiting chebyshev polynomials
- Generalisation to additional detector geometries
