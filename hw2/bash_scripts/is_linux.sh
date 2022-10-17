#!/bin/bash

# Check if subdirectory ”Linux” exists in present directory.
# If yes, print message ”course”. If no, print message ”very
# easy” and create the ”Linux” directory.
# Tip: Google for corresponding flag

if [[ ! -e "Linux" ]]; then
   echo "very easy"
   mkdir "Linux"
else
   echo "course"
fi
