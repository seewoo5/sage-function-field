# sage-function-field

Sage codes for function field related computations.

## Requirements

- SageMath
- ipykernel
- numpy
- scipy

Install them inside Sage shell:

```sh
sage --sh
pip install ipykernel numpy scipy
exit
```

Now you can use these packages with Sage kernel.

## How to use

```sh
sh preparse.sh
```
This will generate preparsed python codes including `__init__.py` where you can import them as
```python
from ff import *
```

You can run test codes under `test` by

```sh
sh run_test.sh
```

## Shanks' bias in function fields

`shanks_bias.ipynb` provides supplementary codes for the examples in the paper [Shanks' bias in function fields](https://arxiv.org/abs/2509.16142).

## Powerful Fibonacci polynomials over finite fields

`fibonacci.ipynb` provides supplementary codes for the examples in the paper [Powerful Fibonacci polynomials over finite fields](https://arxiv.org/abs/2601.02664).

## Ties in function field prime race

`chebyshev_tie.ipynb` provides supplementary codes for the examples in the paper [Ties in function field prime race](https://arxiv.org/abs/2603.21005).