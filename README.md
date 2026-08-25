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
gmp
mpfr
boost
g++ # supporting C++17 or higher
```

On Ubuntu you can install all required packages (except for the compiler) by
running the following commands

```bash
sudo apt update
sudo apt install -y git cmake gnuplot libgmp-dev libmpfr-dev libboost-all-dev
```

## Installation and running

In order to download the repository type in the following command

```bash
git clone https://github.com/JKural/enso-cap.git
cd enso-cap
git switch --detach proof-v3
```

Then, in order to build the programs and run the proof, type

```bash
./proof.sh
```

Note that the script downloads sources, builds them and runs the programs.
This process might take several minutes.
