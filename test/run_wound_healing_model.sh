#!/bin/bash
#
# Script to automate running TestWoundHealingInCrossSectionalGeometry and move the output to this results folder.

# Move to the correct directory
cd ~/BuildChaste/

# Make the test
make -j4 TestWoundHealingInCrossSectionalGeometry

# Run the test
ctest -j8 -R TestWoundHealingInCrossSectionalGeometry$ 

# Move the output to a non-temporary directory
cp -r /tmp/axel/testoutput/WoundHealingModel/CrossSection /media/sf_Chaste/projects/WoundHealingModelling/results/



