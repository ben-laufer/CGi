# GAT_CpG_hyper_hypo.sh
# Ben Laufer

# Takes output from DMRichR and performs statistical testing for CpG annotation enrichments

# You can install GAT through conda (conda install -c bioconda gat) and download bioconda here (https://bioconda.github.io/user/install.html#install-conda)

# The working directory requires:
# 1) A bed file of DMRs (DMRs.bed) from DMRichR (Extra > GAT).
# 2) A bed file of background regions (background.bed) from DMRichR (Extra > GAT).
# 3) A bed file of CpG annotations (hg38_CpG.bed) from DMRichR (Extra > GAT).
# 4) A isochore file that allows you to correct for GC content.
# isochores from: http://bioinfo2.ugr.es:8080/isochores/maps/
# ioschores lifted over from hg19 to hg38 and edited in R for bed format: chr, start, stop, "iso" and then sorted
# Ref:https://github.com/AndreasHeger/gat/issues/4:

# If using external bed files, make sure they have unix and not mac or windows line breaks

echo "Sorting bed files"

sort -k1,1 -k2,2n DMRs.bed > DMRs_sorted.bed
sort -k1,1 -k2,2n background.bed > background_sorted.bed
sort -k1,1 -k2,2n hg38CpG.bed > hg38CpG_sorted.bed

echo "Done sorting"

echo "Running GAT for CpG annotation enrichments"

call="gat-run.py \
--segments=DMRs_sorted.bed \
--with-segment-tracks \
--annotations=hg38CpG_sorted.bed \
--workspace=background_sorted.bed \
--isochore-file=hg38isochores_sorted.bed \
--counter=segment-overlap \
--num-samples=10000 \
--num-threads=10 \
--log=Gat_CpG_hyper_hypo.log \
> Gat_CpG_hyper_hypo_results.tsv"

echo $call
eval $call

echo "Done enriching"

echo "Removing temporary files"
rm DMRs_sorted.bed
rm background_sorted.bed
rm hg38CpG_sorted.bed

echo "Done"
