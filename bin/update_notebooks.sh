#!/bin/bash

for i in $(ls examples/*/*.ipynb ); do
    echo $i; jupyter nbconvert --to notebook $i --output $(basename $i)
done
