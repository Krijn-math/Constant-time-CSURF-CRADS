#!/bin/bash

prime=$1
multieval=$2

( python3 tunned-parameters.py -p $prime -m $multieval -a csidh) &> ijk/$prime-$multieval
