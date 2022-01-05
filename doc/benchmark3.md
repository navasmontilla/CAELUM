# Benchmark #1: 1D Riemann problems

The following Riemann problems are computed using the WENO method in combination with the HLLE solver. The grid is discretized in 100 computational cells. 

These cases can be computed using the Jupyter Notebook [1D_multicomponent.ipynb](../1D_RPs.ipynb), which can executed by running in Anaconda:
``` 
jupyter notebook 1D_multicomponent.ipynb
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
    <td><img src="https://render.githubusercontent.com/render/math?math=\gamma"></td>
    <td>1.4</td>
    <td>1.6</td>
  </tr>
 </table>
 

<figure style="text-align: center;">
  <img src="bench_multicomponent_comp_1.png" width="100%" alt="my alt text"/>
</figure>

