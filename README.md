## Hyperbolic elements in affine Coxeter groups

This program checks Lemma 3.21 of [PS19] for the exceptional affine Coxeter groups F4, E6, E7, E8 (the case G2 can be easily checked by hand).
The lemma is checked for all hyperbolic elements *u* in [1,*w*], without the irreducibility hypothesis.

Requirements: Sage 8.8 and Python 2.7.


### Usage

```bash
sage check_hyperbolic.sage F4
sage check_hyperbolic.sage E6
sage check_hyperbolic.sage E7
sage check_hyperbolic.sage E8
```

The optional argument `-v` can be added to print more information.

The case C is also implemented. For example:

```bash
sage check_hyperbolic.sage C2
```

### Bibliography

[Arm09] D. Armstrong, *Generalized noncrossing partitions and combinatorics of Coxeter groups*, Memoirs of the American Mathematical Society (2009).

[Hum92] J. E. Humphreys, *Reflection groups and Coxeter groups*, Cambridge University Press (1992).

[MS17] J. McCammond and R. Sulway, *Artin groups of Euclidean type*, Inventiones Mathematicae 210(1), 231-282 (2017).

[PS19] G. Paolini and M. Salvetti, *Proof of the K(&#x03C0;,1) conjecture for affine Artin groups*, arXiv preprint (2019).
