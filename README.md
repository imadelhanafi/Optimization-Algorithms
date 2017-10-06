
# Optimization Algorithms

In this project we want to solve numerically an optimization problem with constraints.

 ![equation](http://latex.codecogs.com/gif.latex?%24%24%20%28P_%7BE%7D%29%20%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Brl%7D%20%5Ctext%7Bmin%7D%20%26%20e%28x%2Cy%29%5C%5C%20c_E%28x%2Cy%29%20%26%3D%200%20%5C%5C%20%5Cend%7Barray%7D%20%5Cright.%20%24%24) 

### Structure

This project is devided into two parts : 

* **Optimizer** - Implementation of optimization algorithms : Newton and Generalized Newton using SQP (Sequential Quadratic Programming).
* **Simulator** - Calculate necessary parameters (from a real problem) that will be used as inputs in the optimizer.

### Problem
Our goal is to find the position of static equilibrium of a chain formed of rigid bars in a vertical plane.

Inputs of this problems (defined in ch.m) are:

* LB : Lengths of bars
* xy : position of the nodes
* (a,b) : Position of end node, knowing that the other extremity is at (0,0)
* Max-iteration : max iteration allowed in the optimization algorithms
* Tolerance : Error tolerance allowed
* imode : choice of the optimization algorithms

The following figure represents a chain (black points are nodes):

<p align="center">
  <img src="http://code2net.com/folder/node.png" width="400">
</p>


## Simulator

In this part we will get all the necessary functions and parameters (implemented in **chs.m**) to apply the optimization algorithm:

* e: the potential energy in the position xy
* Ce: the value of constraints (length of bars) in the position xy. The Length of bars should be constant in every step
* g: Gradiant of e at xy
* Ae: Jacobian of contraints at xy
* Hl: the Hessian of the Lagrangian at xy

These parameters are returned by the function *chs.m* using the initila inputs from *ch.m*.


## Optimizer
