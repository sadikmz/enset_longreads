# Mazia

# EV_mazia_trf_1bp.split.bed
# EV_mazia_trf_2bp.split.bed
# EV_mazia_trf_3bp.split.bed
# EV_mazia_trf_4bp.split.bed
# EV_mazia_trf_6bp.split.bed
# EV_mazia_trf_7bp.split.bed
# EV_mazia_trf_8bp.split.bed
# EV_mazia_trf_9bp.split.bed
# EV_mazia_trf_10bp.split.bed
# EV_mazia_trf_17bp.split.bed
# EV_mazia_trf_20bp.split.bed
# EV_mazia_trf_21bp.split.bed
# EV_mazia_trf_25bp.split.bed
# EV_mazia_trf_29bp.split.bed
# EV_mazia_trf_38bp.split.bed
# EV_mazia_trf_14bp.split.bed
# EV_mazia_trf_15p.split.bed
# EV_mazia_trf_19bp.split.bed




EV_mazia_trf_134bp.split.bed
EV_mazia_trf_145bp.split.bed
EV_mazia_trf_146bp.split.bed
EV_mazia_trf_278bp.split.bed
EV_mazia_trf_279bp.split.bed
EV_mazia_trf_280bp.split.bed
EV_mazia_trf_281bp.split.bed
EV_mazia_trf_291bp.split.bed
EV_mazia_trf_292bp.split.bed
EV_mazia_trf_413bp.split.bed
EV_mazia_trf_414bp.split.bed
EV_mazia_trf_423bp.split.bed
EV_mazia_trf_424bp.split.bed
EV_mazia_trf_425bp.split.bed
EV_mazia_trf_426bp.split.bed

#  save the above as mazia_centromere.txt 

for i in $(cat mazia_centromere.txt);
do 
	cat $i >> mazia_centromere.bed
done

bedtools sort -i mazia_centromere.bed | bedtools merge > mazia_centromere.merged.bed 





# EV_wildb