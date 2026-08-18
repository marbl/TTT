#!/usr/bin/python
import pulp

# Check available solvers
print(pulp.listSolvers(onlyAvailable=True))
