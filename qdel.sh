#!/bin/bash

for j in $(seq $1 $2); do
  qdel $j
done

