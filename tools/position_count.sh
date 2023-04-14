#!/bin/bash

# Get the current directory
dir=$(pwd)

# Prompt the user to enter a number starting with SRR or ERR
read -p "Enter the SRR or ERR number: " number

# Search for directories matching the pattern SRR* or ERR*
directories=(${dir}/SRR* ${dir}/ERR*)

# Loop over the directories and find the match
match_found=false
match_position=0
for ((i=0; i<${#directories[@]}; i++))
do
    if [[ ${directories[i]} == ${dir}/${number}* ]]; then
        match_position=$((i+1))
        match_found=true
        break
    fi
done

# Print a message indicating the position of the matching directory
if [ "$match_found" = true ]; then
    echo "Match found at position $match_position out of ${#directories[@]} directories"
else
    echo "No match found for number $number"
fi

