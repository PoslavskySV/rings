#/bin/bash

export algebench_cmd=algebench
export output_path=problems
export nproblems=110

for char in 0 2 524287
do
    for nvars in 3 4 5 6 7
    do
        # generate non-trivial problems
        ${algebench_cmd} generate factor uniform \
                        --n-problems ${nproblems} \
                        --n-variables ${nvars} \
                        --characteristic ${char} \
                        --bit-length 16 \
                        --n-factors 3 \
                        --min-degree 0 \
                        --max-degree 15 \
                        --size 20 \
                        "${output_path}/factor_${nvars}_char_${char}.problems"
        
        # generate trivial problems
        ${algebench_cmd} generate factor uniform \
                        --n-problems ${nproblems} \
                        --n-variables ${nvars} \
                        --characteristic ${char} \
                        --bit-length 16 \
                        --n-factors 3 \
                        --min-degree 0 \
                        --max-degree 15 \
                        --size 20 \
                        --trivial \
                        "${output_path}/factor_${nvars}_char_${char}_trivial.problems"
    done
done
