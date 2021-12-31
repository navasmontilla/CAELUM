# Benchmark #2: 2D Riemann problems

The following Riemann problems are computed using the WENO method in combination with the HLLE solver. The grid is discretized in 160x160 computational cells. 

These cases can be computed using the Jupyter Notebook [2D_RPs.ipynb](../2D_RPs.ipynb), which can executed by running in Anaconda:
``` 
jupyter notebook 2D_RPs.ipynb
```

The initial condition for the Riemann problems consists of piecewise constant data in four equal square regions in which the domain is divided, identified as follows:

<table>
  <tr>
    <td>1</td>
    <td>2</td>
  </tr>
  <tr>
    <td>3</td>
    <td>4</td>
  </tr>
 </table> 

## RP1



<table>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=\rho_1=33.0/62.0,">
<img src="https://render.githubusercontent.com/render/math?math=p_1=0.3,">
<img src="https://render.githubusercontent.com/render/math?math=u_1=4/\sqrt{11},">
<img src="https://render.githubusercontent.com/render/math?math=v_1=0"></td>
    <td><img src="https://render.githubusercontent.com/render/math?math=\rho_2=1.5,">
<img src="https://render.githubusercontent.com/render/math?math=p_2=1.5,">
<img src="https://render.githubusercontent.com/render/math?math=u_2=0,">
<img src="https://render.githubusercontent.com/render/math?math=v_2=0"></td>
  </tr>
    <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=\rho_3=77.0/558.0,">
<img src="https://render.githubusercontent.com/render/math?math=p_3=9.0/310.0,">
<img src="https://render.githubusercontent.com/render/math?math=u_3=4/\sqrt{11},">
<img src="https://render.githubusercontent.com/render/math?math=v_3=4/\sqrt{11}"></td>
    <td><img src="https://render.githubusercontent.com/render/math?math=\rho_4=33.0/62.0,">
<img src="https://render.githubusercontent.com/render/math?math=p_4=0.3,">
<img src="https://render.githubusercontent.com/render/math?math=u_4=4/\sqrt{11},">
<img src="https://render.githubusercontent.com/render/math?math=v_4=0"></td>
  </tr>
 </table> 
 

<figure style="text-align: center;">
  <img src="bench_RP2D_1.png" width="100%" alt="my alt text"/>
</figure>

<figure style="text-align: center;">
  <img src="rp1_art.tif.jpg" width="100%" alt="my alt text"/>
</figure>



## RP2


<table>
  <tr>
    <td></td>
    <td>Left</td>
    <td>Right</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=\rho"></td>
    <td>1.0</td>
    <td>0.125</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=p"></td>
    <td>1.0</td>
    <td>0.1</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=u"></td>
    <td>0.0</td>
    <td>0.0</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=\phi"></td>
    <td>1.666</td>
    <td>5.0</td>
  </tr>
 </table>
 
<figure style="text-align: center;">
  <img src="bench_RP2D_2.png" width="100%" alt="my alt text"/>
</figure>

<figure style="text-align: center;">
  <img src="rp2_art.tif.jpg" width="100%" alt="my alt text"/>
</figure>