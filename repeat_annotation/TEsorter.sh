#!/bin/bash

genotype= #mazia/wildb/wildc/epo
TEsorter \
${genotype}.intact.LTR_TE.fa \
--seq-type nucl \
--processors 48 \
--pass2-rule 80-80-80 \
--prefix ${genotype}_intact.LTR.TE.tesorter.out \
-db rexdb-plant \
-genome
