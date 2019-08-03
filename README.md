# 2d-improved-Fast-Marching-Method
## Solving the Eikonal equation in 2 dimension, with fast marching method.

This repository is a slight modification of the package: scikit-fmm: 
https://github.com/scikit-fmm/scikit-fmm.


The code is changed to use the linear approximation for the situation when the quadratic equation appeared in the update formula has no solution, which is not incorporated by the original scikit-fmm package. 

The 2d numerical approximation formula can be found in:
https://en.wikipedia.org/wiki/Eikonal_equation#Numerical_approximation.

### Requirement
- Python 3.6
- [Numpy](http://www.numpy.org/)
- [Scipy](https://www.scipy.org/)

### Installation
In command line, enter
``` 
python setup.py install
```








