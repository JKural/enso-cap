#!/usr/bin/env bash

proof_executable=build/programs/fixed_point_proof

kappas=(
  "[0.001,0.100]"
  "[0.100,0.150]"
  "[0.150,0.200]"
  "[0.200,0.250]"
  "[0.250,0.270]"
  "[0.270,0.280]"
  "[0.280,0.285]"
  "[0.285,0.290]"
  "[0.290,0.293]"
  "[0.293,0.296]"
  "[0.296,0.298]"
  "[0.298,0.300]"
)

a="[1,1]"
b="[3,3]"
tau="[1.2,1.2]"
p=60
n=2
output_directory="output"

if [[ ! -x "$proof_executable" ]]; then
  echo "Proof executable does not exist. Building..."
  ./build.sh
  echo "Done"
fi

for kappa in "${kappas[@]}"; do
  echo "Proving orbits for kappa = $kappa"
  "$proof_executable" --a="$a" \
    --b="$b" \
    --tau="$tau" \
    --p="$p" \
    --n="$n" \
    --output_directory="$output_directory" \
    --output_name="output-kappa=$kappa" \
    --kappa="$kappa"
  if [[ $? ]]; then
    echo "Proof failed. Aborting..."
    exit 1
  fi
done

echo "Proof successful"
