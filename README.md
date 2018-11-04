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

The example above shows the outputs obtained for the Ishigami function defined earlier and for a sample size s=8000.

### Graphic output

![output](https://user-images.githubusercontent.com/44723660/47966239-7f25f680-e050-11e8-98c5-e341eb33f8fa.jpg)

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

Number of samples        : 8000 
Maximum order of indices : 3 
Number of factors        : 3 
-> x, sampled by LHS according to uniform continuous law of parameters min = -3.141593 and max = 3.141593
-> y, sampled by LHS according to uniform continuous law of parameters min = -3.141593 and max = 3.141593
-> z, sampled by LHS according to uniform continuous law of parameters min = -3.141593 and max = 3.141593

----------- Results of the sensitivity analysis ------------

Sensivity indices of 1st order :
x :  0.311
y :  0.4346
z : -0.03017

Sensivity indices of 2nd order :
x,y :  0.01387
x,z :  0.2706
y,z :  0.01387

Sensivity indice of 3rd order :
x,y,z : -0.01387

Sum of sensitivity indices of 1st order :  0.7155
Sum of sensitivity indices of 2nd order :  0.2984
Sum of sensitivity indices of 3rd order : -0.01387

Calculation by complementarity
Total sensitivity index of x :  0.5817
Total sensitivity index of y :  0.4485
Total sensitivity index of z :  0.2405

Calculation by sum
Total sensitivity index of x :  0.5817
Total sensitivity index of y :  0.4485
Total sensitivity index of z :  0.2405

------------- Cumulative script execution time -------------

Initialization                                  :        0ms
Creation of matrices A and B                    :        0ms
Creation of matrices C                          :     1s26ms
Collection of the outputs Y                     :     1s89ms
Calculation of sensitivity indices              : 43m33s99ms
Estimation of confidence intervals by bootstrap : 43m34s57ms
Creation of the output graph                    : 44m00s20ms
Storage of the info.txt file                    : 44m00s71ms
End                                             : 44m00s71ms
```

### Results comparison 

Indices | Theoriticals | s=1000 | s=2000 | s=4000 | s=8000 | s=16000
------- | ------------ | ------ | ------ | ------ | ------ | -------
S1 | 0.3139 | 0.2822 | 0.2891 | 0.2989 | 0.311 |
S2 | 0.4424 | 0.4341 | 0.4432 | 0.4594 | 0.4346 |
S3 | 0 | -0.03871 | -0.01467 | 0.02521 | -0.03017 |
S12 | 0 | 0.01615 | 0.01901 | 0.0019 | 0.01387 |
S13 | 0.2437 | 0.3063 | 0.2634 | 0.2146 | 0.2706 |
S23 | 0 | 0.01615 | 0.01901 | 0.0019 | 0.01387 |
S123| 0 | -0.01615 | -0.01901 | -0.0019 | -0.01387 |
ST1 | 0.5576 | 0.5885 | 0.5525 | 0.5135 | 0.5817 |
ST2 | 0.4424 | 0.4502 | 0.4622 | 0.4613 | 0.4485 |
ST3 | 0.2437 | 0.2676 | 0.2487 | 0.2398 | 0.2405 |

Distances from theoriticals results | s=1000 | s=2000 | s=4000 | s=8000 | s=16000
----------------------------------- | ------ | ------ | ------ | ------ | -------
Indices of x-th order | 0.0853 | 0.0480 | 0.0448 | 0.0478 | 
Total order indices | 0.0398 | 0.0210 | 0.0481 | 0.0251 |
