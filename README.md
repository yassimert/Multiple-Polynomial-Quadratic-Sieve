# Multiple-Polynomial-Quadratic-Sieve
C Implementation of the Multiple Polynomial Quadratic Sieve Algorithm

## How to Compile & Run

```
gcc -O3 main.c mpqs.c -o mpqs.out -lgmp -lm && ./mpqs.out
```

## Requirements
- [GNU Multiple Precision Arithmetic Library](https://gmplib.org)
- A local [Magma Computer Algebra](http://magma.maths.usyd.edu.au/magma/) installation (to be able to factorize larger numbers) 
- [Python3](https://www.python.org/downloads/) for _encoder.py_. If you have a local Magma installed, there is no need for _encoder.py_
