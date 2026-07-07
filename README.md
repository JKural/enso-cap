# A computer assisted proof of existence of periodic orbit in ENSO DDE model

This repository implements a computer assisted proof of existence of periodic orbits
in a DDE model of ENSO given by

$$
  h'(t) = -a \tanh(\kappa h(t-\tau)) + b \cos(2 \pi t)
$$

## Dependencies

```bash
git
cmake
gnuplot
g++ # supporting C++17 or higher
```

## Installation and running

In order to download the repository type is the following command

```bash
git clone https://github.com/JKural/enso-cap.git
cd enso-cap
git switch --detach proof-v2
```

Then, in order to build the programs and run the proof, type

```bash
./proof.sh
```
