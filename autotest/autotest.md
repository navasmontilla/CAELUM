## AUTOMATIC TEST RESULTS (autotest.py)

### Convergence rates

The linear transport equation is solved in 1D inside a cyclic domain x=[0,1] using the following initial conditions:

$$u(x,0)=1+0.5\sin(2\pi x)$$

The advection velocity is $a=1$ and the simulation time is t=5 s. The numerical errors are computed using the L_1 error norms by comparing with the exact solution:

$$u(x,0)=1+0.5\sin(2\pi (x-at))$$

<figure style="text-align: center;">
  <img src="caseLinear/out/convergences.png" width="50%" alt="my alt text"/>
</figure>
