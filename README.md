# Sensivity Analysis

This script allows you to perform a sensitivity analysis for a given function, with Sobol method. The parameters of the function to analyze are sampled according to the LHS method, considering that each parameter is derived from a continuous uniform law. It is therefore necessary to know the range of variation of each parameter. Each of the calculated indices is framed using the bootstrap method.

## Running the tests

The given script is configured for the example of the Ishigami function with a=7 and b=0.1.
The sensivity indices obtained approximate the theoretical results, namely:
```
S1  : 0.3139
S2  : 0.4424
S3  : 0
S12 : 0
S23 : 0
S13 : 0.2437
S123: 0
ST1 : 0.5576
ST2 : 0.4424
ST3 : 0.2437
```

## Outputs

The script manages the post-processing of sensitivity indices by creating a graphic that represents first-order sensitivity indices.
A raw text file is also created and summarizes the script inputs, x-th order sensitivity indices, total sensitivity indices as well as the intermediate execution times of the script steps.

### Graphic output

![output](https://user-images.githubusercontent.com/44723660/47964762-98be4280-e03e-11e8-8c1a-ecab8957d8aa.jpg)

### Raw output

```
------------------------------------------------------------
------------------- Sensitivity analysis -------------------
------------------------------------------------------------

------------------- Function parameters --------------------

function(x,y,z){
  return(sin(x)+7*sin(y)^2+0.1*z^4*sin(x))
}

------------- Sensitivity analysis parameters --------------

Number of samples        : 7000 
Maximum order of indices : 3 
Number of factors        : 3 
-> x, sampled by LHS according to uniform continuous law of parameters min = -3.141593 and max = 3.141593
-> y, sampled by LHS according to uniform continuous law of parameters min = -3.141593 and max = 3.141593
-> z, sampled by LHS according to uniform continuous law of parameters min = -3.141593 and max = 3.141593

----------- Results of the sensitivity analysis ------------

Sensivity indices of 1st order :
x :  0.3119
y :  0.4253
z : -0.002792

Sensivity indices of 2nd order :
x,y :  0.002839
x,z :  0.2627
y,z :  0.002839

Sensivity indice of 3rd order :
x,y,z : -0.002839

Sum of sensitivity indices of 1st order :  0.7344
Sum of sensitivity indices of 2nd order :  0.2684
Sum of sensitivity indices of 3rd order : -0.002839

Calculation by complementarity
Total sensitivity index of x :  0.5746
Total sensitivity index of y :  0.4281
Total sensitivity index of z :  0.26

Calculation by sum
Total sensitivity index of x :  0.5746
Total sensitivity index of y :  0.4281
Total sensitivity index of z :  0.26

------------- Cumulative script execution time -------------

Initialization                                  :        0ms
Creation of matrices A and B                    :        0ms
Creation of matrices C                          :     1s06ms
Collection of the outputs Y                     :     1s82ms
Calculation of sensitivity indices              : 31m57s24ms
Estimation of confidence intervals by bootstrap : 31m57s65ms
Creation of the output graph                    : 32m20s33ms
Storage of the info.txt file                    : 32m20s64ms
End                                             : 32m20s65ms
```
