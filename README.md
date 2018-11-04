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

The script manages the post-processing of sensitivity indices by creating a visual that represents first-order sensitivity indices.
A text file is also created and summarizes the script inputs, x-th order sensitivity indices, total sensitivity indices as well as the intermediate execution times of the script steps.
