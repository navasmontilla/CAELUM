# Benchmark #4: Convergence rate test

The Euler equations are computed inside a cyclic domain x=[0,1] using the following initial conditions:

<img src="https://latex.codecogs.com/svg.image?\left\{&space;\begin{array}{l}p=1&space;\\u=1&space;\\u=w=0&space;\\\rho=1&plus;0.5\sin(2\pi&space;x)\end{array}\right.&space; "/>

setting up a linear transport problem. The simulation time is t=5 s. The numerical errors are computed using the L_1 and L_inf error norms by comparing with the exact solution:

<img src="https://latex.codecogs.com/svg.image?\left\{&space;\begin{array}{l}p=1&space;\\u=1&space;\\u=w=0&space;\\\rho=1&plus;0.5\sin(2\pi&space;(x-ut))\end{array}\right.&space; "/>

### h-refinement

<figure style="text-align: center;">
  <img src="convergence_href.png" width="100%" alt="my alt text"/>
</figure>

### p-refinement

<figure style="text-align: center;">
   <img src="convergence_pref.png" width="100%" alt="my alt text"/>
</figure>




