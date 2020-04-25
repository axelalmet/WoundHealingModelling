#!/bin/bash


# Move to the correct directory
cd ~/BuildChaste/

# Make the test
make $1

# Run the test
ctest -R $1

# Move the output to a non-temporary directory
cp -r /tmp/axel/testoutput/WoundHealingModel/CrossSection $2

# Move LastTestLog to a non-temporary directory (this is useful for debugging)
cp -r ~/BuildChaste/Testing/Temporary/LastTest.log $2