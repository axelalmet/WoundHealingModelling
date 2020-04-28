#!/bin/bash

# Move to the correct directory
cd ~/BuildChaste/

# Make the test
CUDA_VISIBLE_DEVICES=$1 make $2

# Run the test
CUDA_VISIBLE_DEVICES=$1 ctest -R $2

# Move the output to a non-temporary directory
cp -r /tmp/axela/testoutput/WoundHealingModel/CrossSection $3

# Move LastTestLog to a non-temporary directory (this is useful for debugging)
cp -r ~/BuildChaste/Testing/Temporary/LastTest.log $3