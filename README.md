
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
* Tolerance : Error tolerance
* imode : choice of the optimization algorithms

The following figure represents a chain (black points are nodes):

<p align="center">
  <img src="http://code2net.com/folder/screenshots/node.png" width="400">
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

Optimization algorithms are implemented in *sqp.m*

### Newton Method

Let *f* be the function that we want to minimize and *c* the constraints.

 ![equation](http://latex.codecogs.com/gif.latex?%24%24%20%28P%29%20%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Brl%7D%20%5Ctext%7Bmin%7D%20%26%20f%28x%29%5C%5C%20c%28x%29%20%26%3D%200%20%5C%5C%20%5Cend%7Barray%7D%20%5Cright.%20%24%24) 
 
 We want to find a primal-dual solution ![equation](http://latex.codecogs.com/gif.latex?%24%24%20%28x_*%2C%5Clambda_*%29%20%24%24) that satisfy the following optimality condition:
 
 ![equation](http://latex.codecogs.com/gif.latex?%24%24%20%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Brl%7D%20%5Cnabla%20f%28x_*%29&plus;c%27%28x_*%29%5ET%5Clambda_*%20%26%20%3D%200%5C%5C%20c%28x_*%29%20%26%3D%200%20%5C%5C%20%5Cend%7Barray%7D%20%5Cright.%20%24%24)
 
Solving this problem is equivalent to find the zero of the following function by  setting ![equation](http://latex.codecogs.com/gif.latex?%24%24%20z%20%3D%20%28x%2C%5Clambda%29%20%24%24) :

![equation](http://latex.codecogs.com/gif.latex?%24%24%20F%28z%29%3D%5Cbegin%7Bpmatrix%7D%20%5Cnabla%20f%28x%29&plus;c%27%28x%29%5ET%5Clambda%20%5C%5C%20c%28x%29%20%5Cend%7Bpmatrix%7D%20%24%24)

To find zero z* of this function we will use Newton Algorithm and at iteration z_k we need to solve the following linear system (finding x_k and lambda_k)

![equation](http://latex.codecogs.com/gif.latex?%24%24%20%5Cbegin%7Barray%7D%7Brl%7D%20%5Cbegin%7Bpmatrix%7D%20L_k%20%26%20A_k%5ET%20%5C%5C%20A_k%20%26%200%20%5C%5C%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20d_k%5C%5C%20%5Clambda_k%5E%7BPQ%7D%20%5C%5C%20%5Cend%7Bpmatrix%7D%20%26%3D%20-%20%5Cbegin%7Bpmatrix%7D%20%5Cnabla%20f%28x_k%29%5C%5C%20c%28x_k%29%20%5Cend%7Bpmatrix%7D%20%5Cend%7Barray%7D%20%24%24) 
with ![equation](http://latex.codecogs.com/gif.latex?%24L_k%3D%5Cnabla_%7Bxx%7Dl%28x_k%2C%5Clambda_k%29%24) the Hessian of the Lagrangian and ![equation](http://latex.codecogs.com/gif.latex?%24A_k%3Dc%27%28x_k%29%24 ).

Then:

![equation](http://latex.codecogs.com/gif.latex?%5C%5C%20%24x_%7Bk&plus;1%7D%3Dx_k&plus;d_k%24%24%5C%5C%20%24%5Clambda_%7Bk&plus;1%7D%3D%5Clambda_k%5E%7BPQ%7D%24)
 
**Convergance** -
**Nature of the solution**

### Generalized Newton Method (Sequential Quadratic Programming)


The aim of this part is to modify the newton method so that it converges to a solution whatever the initial position.

To solve this problem we will use a penalty function.

Penalization is a concept that transforms an optimization problem with constraints into a sequence of optimization problems without constraint. The goal is to introduce a new function whose value is penalized when the constraints are not respected.


**Principle**

In Newton's algorithm we have that:
![equation](http://latex.codecogs.com/gif.latex?%5C%5C%20%24x_%7Bk&plus;1%7D%3Dx_k&plus;d_k%24%24%5C%5C%20%24%5Clambda_%7Bk&plus;1%7D%3D%5Clambda_k%5E%7BPQ%7D%24)

Such that  ![equation](http://latex.codecogs.com/gif.latex?%24%24%28d_k%2C%5Clambda_k%5E%7BPQ%7D%29%24%24) is primal-dual solution of the following quadratic osculator problem ([see SQP](https://en.wikipedia.org/wiki/Sequential_quadratic_programming))

![equation](http://latex.codecogs.com/gif.latex?%24%24%20%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Brl%7D%20%5Ctext%7Bmin%7D_%7Bd%7D%20%5Cnabla%20f%28x_k%29%5ETd%20&plus;%20%5Cfrac%7B1%7D%7B2%7Dd%5ETH_kd%20%26%5C%5C%20c%28x_k%29&plus;c%27%28x_k%29d%20%26%3D%200%20%5C%5C%20%5Cend%7Barray%7D%20%5Cright.%20%24%24)

Hk is the hessian of the lagrangian. To solve this problem in a reasonable time, the value of the hessian is modified slightly so that it becomes a positive matrix. We thus use the Cholesky factorization, which makes it possible to find a matrix Mk defined positive close to Hk.

We introduce the penalization function:

![equation](http://latex.codecogs.com/gif.latex?%24%5Ctheta_%7B%5Csigma%7D%28x%29%3Df%28x%29&plus;%5Csigma%20%5C%7C%20c%28x%29%20%5C%7C%20_1%24)

We know that:
![equation](http://latex.codecogs.com/gif.latex?%24%5Ctheta%27_%7B%5Csigma%7D%28x_k%2Cd_k%29%3D-d_k%5ETM_kd_k&plus;%28%5Clambda_k%5E%7BPQ%7D%29%5ETc_k-%5Csigma%20%5C%7C%20c%28x%29%20%5C%7C%20_1%24)

Line search will be made on ![equation](http://latex.codecogs.com/gif.latex?%24%5Ctheta_%7B%5Csigma%7D%24)  and the goal now is to find alpha such that

![equation](http://latex.codecogs.com/gif.latex?%24x_%7Bk&plus;1%7D%20%3Dx_%7Bk%7D&plus;%5Calpha_%7Bk%7Dd_%7Bk%7D%24) and
![equation](http://latex.codecogs.com/gif.latex?%24%24%20%5Ctheta_%7B%5Csigma%7D%28x_k&plus;%5Calpha_k%20x_k%29%20%5Cleq%20%5Ctheta_%7B%5Csigma%7D%28x_k%29%20&plus;%20%5Comega%20%5Calpha_k%20%5CDelta_k%20%24%24)  
with ![equation](http://latex.codecogs.com/gif.latex?%24%5CDelta_k%3D-d_k%5ETM_kd_k%20&plus;%20%28%5Clambda_k%5E%7BPQ%7D%29%5ETc_k-%5Csigma%20%5C%7C%20c%28x%29%20%5C%7C%20_1%24.%5C%5C) and ω=0.0001.

Cholesky decomposition is in: 
Implementation of this method is in:



