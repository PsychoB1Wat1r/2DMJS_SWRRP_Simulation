## Networked Control of 2DMJS via SWRRP

This repository contains the code implementation of our patent (Protocol-based AETC for heat exchangers).

## Explanation:
We use **MATLAB** as the programming language and employ the **Mosek solver** along with the **YALMIP toolbox** as development tools to perform numerical solution and simulation for the proposed method.

## Mosek Solver:
Mosek is a high-performance mathematical optimization solver specifically designed for convex optimization problems. For official documentation, please refer to this [link](https://www.mosek.com/documentation/).

## YALMIP Toolbox:
YALMIP (Yet Another LMI Parser) is an optimization modeling toolbox for MATLAB that facilitates the formulation and solution of various optimization problems. Developed by [**Johan Löfberg**](https://scholar.google.com/citations?user=No-9sDUAAAAJ&hl=en), YALMIP is primarily used for convex optimization but can also handle certain non-convex problems.

It is important to note that YALMIP itself is **not a solver**; rather, it serves as a modeling interface that translates user-defined optimization problems into a standard form and then calls external solvers such as **Mosek**, **Gurobi**, **SDPT3**, **SeDuMi**, and others for solution. For official documentation, please refer to this [link](https://yalmip.github.io/).
