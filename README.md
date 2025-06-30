# Factoring through Frobenius

Code accompaining the paper "Orient Express: Using Frobenius to Express
Oriented Isogenies".

## Square-free case

To test key exchange in the square-free case run

```
sage --python -O cs_ftf.py
```

Different security levels can be selected inside the file. To run timings for
all security level run

```
sage --python -O cs_timings.py
```

To generate new parameters, go to `/params` and run
```
sage cs_pgen.py
```

## Square case

To test key exchange in the square case run

```
sage --python -O sub_ftf.py
```

Different security levels can be selected inside the file. To run timings for
all security level run

```
sage --python -O sub_timings.py
```

To generate new parameters, go to `/params` and run
```
sage --python sub_pgen.py B a f1 f2 m n n_primes
```
