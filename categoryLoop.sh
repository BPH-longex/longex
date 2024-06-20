#!/bin/bash

for i in {0..7}
do
    root '$1($i)' &
done