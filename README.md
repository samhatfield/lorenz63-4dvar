# lorenz63-4dvar
The 4DVar data assimilation technique implemented with the Lorenz '63 model

This repository contains an implementation of the 4DVar technique with the Lorenz '63 model. The code is intended to be a simple introduction to this technique, and includes both the "plain" formulation of the 4DVar technique, as well as the incremental form (in the "incremental" branch).

Minimisation is currently performed with a very simple gradient descent method. In the future I would like to implement a conjugate gradient technique.

## Install
Simply `git clone` the repository. Then run `make` to build the system and `./main` to run it. The output will be stored in a serious of plaintext files. The output can be visualised by running `python plot.py`. This will require matplotlib and numpy.

## Tests
Some tests for the tangent linear and adjoint models are included in the `test` repository. To build these, first make sure to run `make clean` in the parent directory (the one containing `test`). Then `cd` to `test` and run `make` to build and `./test` to run.
