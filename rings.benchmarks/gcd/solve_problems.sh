#!/bin/bash

export algebench_cmd=algebench
export input_path="problems"
export output_path="results"

export singular_exec="/Applications/Singular.app/Contents/bin/Singular"
export mathematica_exec="/Applications/Mathematica.app/Contents/MacOS/wolframscript"
export fermat_exec="/Users/poslavskysv/Downloads/ferm6/fer64"
export form_exec="form"

function solve(){
    input=$1
    output="${input%.*}"
    output="$output.results.tsv"
    ${algebench_cmd} solve \
            --threads 2 \
            --rings \
            --mathematica --mathematica-exec ${mathematica_exec} \
            --singular --singular-exec ${singular_exec} \
            --form --form-exec ${form_exec} \
            --fermat --fermat-exec ${fermat_exec} \
            "${input_path}/$input" "${output_path}/$output"
}

export -f solve

ls -1 ${input_path} | grep ".problems$" | parallel -j1 'solve {}'