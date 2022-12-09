#!/bin/bash

Files=(
    thesis.tex \
    intro.tex \
    design/aistats/bayesian.tex \
    design/continuous_nips/continuous.tex \
    design/equivalents_nips/neurips_2020.tex \
    ic/aistats2021/paper.tex \
    ftvi/ftvi_icml/ftvi_icml.tex \
    hta/neurips_2022.tex
)
for f in ${Files[@]}; do
    aspell \
        --lang=en \
        --mode=tex \
        --conf=./tex_conf \
        --add-extra-dicts=./aspelling_dict \
        check \
        $f
done;