# Scalability and Communication Overhead in Distributed N-Body Simulation

This repository contains the design, implementation, and performance evaluation of a parallel **N-Body simulation** using the **Message Passing Interface (MPI)**. The project explores the scalability of the direct method $O(N^2)$ across various **Google Cloud Platform (GCP)** configurations.

## Overview
The simulation models a closed system of $N$ point masses interacting solely through gravitational forces. To ensure physical stability over long durations, the project utilizes the **semi-implicit Euler integration method**, which is a symplectic integrator that bounds energy oscillations.

### Key Objectives
* Analyze a serial algorithm for solving the N-body problem.
* Perform an a priori study of available parallelism using Amdahl's Law.
* Develop a parallel implementation using **MPI**.
* Evaluate performance and scalability through experiments on **GCP**.

---

## Technical Details

### Physical Model & Complexity
* **Force Calculation**: Based on Newton’s law of universal gravitation and his second law of motion.
* **Time Complexity**: $O(N^2)$ for the direct method, dominated by the acceleration computation loop.
* **Space Complexity**: $O(N)$, requiring storage for masses and 3D vectors for position, velocity, and acceleration.

### Parallel Strategy
The implementation distributes $N$ bodies across $P$ processes.
* **Communication**: Employs `MPI_Allgatherv` to ensure every process has a complete global view of particle positions required for force computation at each step.
* **Custom Types**: Implements a custom MPI data type (`MPI_VECTOR`) to wrap C structs and handle 3D vector data efficiently.
* **Reproducibility**: Uses global index-based seeding for the random number generator to ensure parallel results match serial baselines.

---

## Performance Evaluation
Testing was conducted on multiple GCP configurations, from single 16-vCPU "fat" instances to clusters distributed across four geographical regions.

### Scaling Results
* **Strong Scalability**: Achieved near-ideal linear speedup when processes remain within a single machine or intra-regional clusters where communication overhead is negligible.
* **Inter-regional Bottlenecks**: Performance deteriorates significantly in multi-region deployments because network latency and the blocking nature of `MPI_Allgatherv` dominate the runtime.
* **Weak Scalability**: The simulation exhibits poor weak scalability because the computational workload per process grows with $P$, even if the local data workload is constant.

---

## Repository Structure
* `serial.c`: The baseline serial implementation in C.
* `parallel.c`: The MPI-based parallel implementation.
* `tester.c`: A validation utility to compare serial and parallel output files.
* `serial_testing.c`: A variant for qualitative testing with custom configuration files.
* `animate-nbody_2d.py`: Python script used to generate 2D animations of trajectories.
