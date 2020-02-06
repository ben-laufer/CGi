# GAT.sh
# Perform enrichment testing for CpG and genic annotations within DMRs
# Ben Laufer

echo "Sorting DMR and background region bed files"
sort -k1,1 -k2,2n DMRs.bed > DMRs_sorted.bed
sort -k1,1 -k2,2n background.bed > background_sorted.bed
echo "Done sorting DMR and background region bed files"

#######
# CpG #
#######
echo "Testing for CpG annotation enrichments"

echo "Sorting CpG annotation bed file"
sort -k1,1 -k2,2n hg38CpG.bed > hg38CpG_sorted.bed
echo "Done sorting CpG annotation bed file"

echo "Running GAT for CpG annotation enrichments"
call="gat-run.py \
--segments=DMRs_sorted.bed \
--with-segment-tracks \
--annotations=hg38CpG_sorted.bed \
--workspace=background_sorted.bed \
--isochore-file=hg38isochores.bed \
--counter=segment-overlap \
--num-samples=10000 \
--num-threads=10 \
--log=Gat_CpG_hyper_hypo.log \
> Gat_CpG_hyper_hypo_results.tsv"

echo $call
eval $call
echo "Done enriching for CpG annotations"

echo "CpG testing complete"

#########
# Genic #
#########
echo "Testing for genic annotation enrichments"

echo "Sorting, removing invalid records, and merging genic annotation bed files"
sort -k1,1 -k2,2n enhancers.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > enhancers_merged.bed
sort -k1,1 -k2,2n promoters.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > promoters_merged.bed
sort -k1,1 -k2,2n introns.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > introns_merged.bed
sort -k1,1 -k2,2n boundaries.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > boundaries_merged.bed
sort -k1,1 -k2,2n intergenic.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > intergenic_merged.bed
sort -k1,1 -k2,2n exons.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > exons_merged.bed
sort -k1,1 -k2,2n fiveUTRs.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > fiveUTRs_merged.bed
sort -k1,1 -k2,2n threeUTRs.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > threeUTRs_merged.bed
sort -k1,1 -k2,2n onetofivekb.bed | awk -v OFS='\t' '{if ($2 >= 0) print $1,$2,$3}' | bedtools merge -i stdin > onetofivekb_merged.bed
echo "Done sorting, removing invalid records, and merging bed files"

echo "Running GAT for genic annotation enrichments"
call="gat-run.py \
--segments=DMRs_sorted.bed \
--with-segment-tracks \
--annotations=enhancers_merged.bed \
--annotations=promoters_merged.bed \
--annotations=introns_merged.bed \
--annotations=boundaries_merged.bed \
--annotations=intergenic_merged.bed \
--annotations=exons_merged.bed \
--annotations=fiveUTRs_merged.bed \
--annotations=threeUTRs_merged.bed \
--annotations=onetofivekb_merged.bed \
--workspace=background_sorted.bed \
--isochore-file=hg38isochores.bed \
--counter=segment-overlap \
--num-samples=10000 \
--num-threads=10 \
--log=Gat_genic_hyper_hypo.log \
> Gat_genic_hyper_hypo_results.tsv"

echo $call
eval $call
echo "Done enriching for genic annotations"

echo "Removing temporary files"
rm *_sorted.bed
rm *_merged.bed
echo "Done removing temporary files"

#########
# Plots #
#########

echo "Plotting enrichment results"
call="Rscript \
--vanilla \
GATplots.R"

echo $call
eval $call
echo "Done plotting enrichment results"

echo "Script complete"
