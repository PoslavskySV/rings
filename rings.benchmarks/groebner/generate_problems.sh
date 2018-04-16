#/bin/bash


export algebench_cmd=algebench
export output_path=problems
export KATSURA_SEQ=`seq 5 12`
export CYCLIC_SEQ=`seq 5 9`

for char in 0 1000003
do
    for k in ${KATSURA_SEQ}
    do
        ${algebench_cmd} generate groebner \
                            --problems "katsura_${k}" \
                            --characteristic ${char} \
                            --order grevlex \
                            "${output_path}/katsura_${k}_char_${char}_grevlex.problems"
    done

    for k in ${CYCLIC_SEQ}
    do
        ${algebench_cmd} generate groebner \
                            --problems "cyclic_${k}" \
                            --characteristic ${char} \
                            --order grevlex \
                            "${output_path}/cyclic_${k}_char_${char}_grevlex.problems"        
    done
done
