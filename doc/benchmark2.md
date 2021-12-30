# Benchmark #2: 2D Riemann problems

The following Riemann problems are computed using the WENO method in combination with the HLLE solver. The grid is discretized in 100 computational cells. 

These cases can be computed using the Jupyter Notebook [2D_RPs.ipynb](../2D_RPs.ipynb), which can executed by running in Anaconda:
``` 
jupyter notebook 1D_RPs.ipynb
```

## RP1

<table>
  <tr>
    <td></td>
    <td>Left</td>
    <td>Right</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=\rho"></td>
    <td>1.0</td>
    <td>1.0</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=p"></td>
    <td>1.0</td>
    <td>1.0</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=u"></td>
    <td>0.0</td>
    <td>0.0</td>
  </tr>
  <tr>
    <td><img src="https://render.githubusercontent.com/render/math?math=\phi"></td>
    <td>1.0</td>
    <td>0.0</td>
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