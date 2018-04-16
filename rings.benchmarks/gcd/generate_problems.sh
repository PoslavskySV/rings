#!/bin/bash

export algebench_cmd=algebench
export output_path=problems

### Generating common problems
nproblems=110
for nvars in 3 4 5 6 7 8
do

    for modulus in 0 2
    do

        echo "Generating uniform gcd problems #variables = ${nvars} characteristic = $modulus"

        $algebench_cmd generate gcd uniform \
            --n-problems ${nproblems} \
            --n-variables $nvars \
            --characteristic $modulus \
            --bit-length 32 \
            --size 40 \
            --min-degree 0 --max-degree 30 \
            "${output_path}/gcd_uniform_nvars_${nvars}_characteristic_${modulus}.problems"

        # generate sharp problems only for char 0
        if [ "$modulus" == "0" ]; then
            echo "Generating sharp   gcd problems #variables = $nvars characteristic = $modulus"

            $algebench_cmd generate gcd sharp \
                --n-problems ${nproblems} \
                --n-variables $nvars \
                --characteristic $modulus \
                --bit-length 32 \
                --size 40 \
                --total-degree 50 \
                "${output_path}/gcd_sharp_nvars_${nvars}_characteristic_${modulus}.problems"
        fi

    done
done


### Generating huge problems #vars == 3
for size in 50 100 500 1000 5000 10000
do
    nprobs=10
    if [ "$size" -gt "999" ]; then
        nprobs=3
    fi

    echo "Generating sharp huge gcd problems size = $size (#probs = $nprobs) for #vars = 3"
    
    $algebench_cmd generate gcd sharp \
        --n-problems $nprobs \
        --n-variables 3 \
        --characteristic 0 \
        --bit-length 256 \
        --size $size \
        --total-degree 50 \
        "${output_path}/gcd_huge_${size}_sharp_nvars_3_characteristic_0.problems"
done


### Generating huge problems #vars == 4
for size in 50 100 500 1000 5000
do
    nprobs=3
    if [ "$size" -gt "999" ]; then
        nprobs=1
    fi

    echo "Generating sharp huge gcd problems size = $size (#probs = $nprobs) for #vars = 4"
    
    $algebench_cmd -Xmx8g generate gcd sharp \
        --n-problems $nprobs \
        --n-variables 4 \
        --characteristic 0 \
        --bit-length 256 \
        --size $size \
        --total-degree 50 \
        "${output_path}/gcd_huge_${size}_sharp_nvars_4_characteristic_0.problems"
done
