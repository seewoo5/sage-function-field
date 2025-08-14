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

## Shanks bias in function fields

`shanks_bias.ipynb` provides supplementary codes for the examples in the paper "Shanks bias in function fields".