# Benchmark #4: Convergence rate test

The following Riemann problems are computed using the WENO method in combination with the HLLE solver. The grid is discretized in 100 computational cells. The multi-fluid solver is activated with the following macros:

https://latex.codecogs.com/svg.image?\left\{&space;\begin{array}{l}p=1&space;\\u=1&space;\\u=w=0&space;\\\rho=1&plus;0.5\sin(2\pi&space;x)\end{array}\right.&space;

## h-refinement

<figure style="text-align: center;">
  <img src="convergence_href.png" width="100%" alt="my alt text"/>
</figure>

## p-refinement

<figure style="text-align: center;">
   <img src="convergence_pref.png" width="100%" alt="my alt text"/>
</figure>




