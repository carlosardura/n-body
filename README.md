# N-Body Problem Simulator

This repository contains a three-dimensional simulation of the gravitational interaction of an N-body system in Python.

The generation of the initial conditions follows the **Plummer spherical density model**, designed for the simulation of star clusters in virial equilibrium.

The simulator implements a **Barnes-Hut algorithm** to group distant bodies, approximating their gravitational pull as a single pseudo-particle located at the group’s center of mass, reducing the computational complexity of the force calculation from $O(N^2)$ to $O(NlogN)$.

The algorithm recursively divides the set of bodies into groups by storing them in an octree structure, where each node represents a smaller region of the three-dimensional space. Interactions are evaluated through a **Multipole Acceptance Criterion** (MAC). Collapses due to particle overlap during the tree construction are avoided through a defined maximum recursion depth.

The time evolution of the system is driven by a **Velocity Verlet integrator**. This symplectic numerical method updates, at discrete time steps, the kinematic state of the defined bodies, ensuring long-term energy conservation and stability for Hamiltonian gravitational systems.