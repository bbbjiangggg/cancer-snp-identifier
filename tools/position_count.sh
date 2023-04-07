#!/bin/bash

# Get the current directory
dir=$(pwd)

# Prompt the user to enter a number starting with SRR or ERR
read -p "Enter the number: " number

# Search for directories matching the pattern SRR* or ERR*
directories=(${dir}/SRR* ${dir}/ERR*)

# Loop over the directories and find the match
match_found=false
for ((i=0; i<${#directories[@]}; i++))
do
    if [[ ${directories[i]} == ${dir}/${number}* ]]; then
        echo "Match found at position $i"
        match_found=true
        break
    fi
done

# Print a message if no match was found
if [ "$match_found" = false ]; then
    echo "No match found for number $number"
fi

