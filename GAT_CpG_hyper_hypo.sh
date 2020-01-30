# GAT_CpG_hyper_hypo.sh

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
