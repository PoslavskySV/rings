#!/bin/bash

export algebench_cmd=algebench
export input_path="problems"
export output_path="results"

export singular_exec="/Applications/Singular.app/Contents/bin/Singular"
export mathematica_exec="/Applications/Mathematica.app/Contents/MacOS/wolframscript"

function solve(){
    input=$1
    output="${input%.*}"
    output="$output.results.tsv"
    ${algebench_cmd} solve \
            --rings \
            --singular --singular-exec ${singular_exec} \
            --mathematica --mathematica-exec ${mathematica_exec} \
            "${input_path}/$input" "${output_path}/$output"
}

export -f solve

ls -1 ${input_path} | grep ".problems$" | parallel -j1 'solve {}'