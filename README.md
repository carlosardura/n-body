# N-Body Problem Simulator

This repository contains a three-dimensional simulation of the gravitational interaction of an N-body system in a hybrid Python and C++ architecture.

The generation of the initial conditions follows the **Plummer spherical density model**, designed for the simulation of star clusters in virial equilibrium.

The simulator implements a **Barnes-Hut algorithm** to group distant bodies, approximating their gravitational pull as a single pseudo-particle located at the group’s center of mass, reducing the computational complexity of the force calculation from $O(N^2)$ to $O(NlogN)$.

The algorithm recursively divides the set of bodies into groups by storing them in an octree structure, where each node represents a smaller region of the three-dimensional space. Interactions are evaluated through a **Multipole Acceptance Criterion** (MAC). Collapses due to particle overlap during the tree construction are avoided through a defined maximum recursion depth.

The time evolution of the system is driven by a **Velocity Verlet integrator**. This symplectic numerical method updates, at discrete time steps, the kinematic state of the defined bodies, ensuring long-term energy conservation and stability for Hamiltonian gravitational systems.

To maximize efficiency, the system decouples data orchestration from high-performance computation through an offline execution pipeline.

A **Python frontend** driver defines the global static constants used in the simulation. It manages initialization, generating an **HDF5 file** containing the initial particle conditions defined by the Plummer model. The subsequent execution is then delegated to a C++ engine via a blocking system call, ensuring complete memory isolation between the two environments.

The **C++ backend**, compiled and managed via **CMake** to handle hybrid architecture and HDF5 dependencies, reads the initialization file in order to perform spatial partitioning and numerical integration and stream the output data back to disk. Finally, control is returned to Python upon calculation completion, enabling further data analysis and visualization.