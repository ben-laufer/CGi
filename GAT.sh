# GAT.sh
# Perform enrichment testing for CpG and genic annotations within DMRs
# Ben Laufer

#######
# CpG #
#######
echo "Testing for CpG annotation enrichments"

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
--isochore-file=hg38isochores.bed \
--counter=segment-overlap \
--num-samples=10000 \
--num-threads=10 \
--log=Gat_CpG_hyper_hypo.log \
> Gat_CpG_hyper_hypo_results.tsv"

echo $call
eval $call
echo "Done enriching"

echo "Removing temporary files"
rm *_sorted.bed
echo "Done removing temporary files"

echo "CpG testing complete"

#########
# Genic #
#########
echo "Testing for genic annotation enrichments"

echo "Sorting bed files"
sort -k1,1 -k2,2n DMRs.bed > DMRs_sorted.bed
sort -k1,1 -k2,2n background.bed > background_sorted.bed
sort -k1,1 -k2,2n enhancers.bed > enhancers_sorted.bed
sort -k1,1 -k2,2n promoters.bed > promoters_sorted.bed
sort -k1,1 -k2,2n introns.bed > introns_sorted.bed
sort -k1,1 -k2,2n boundaries.bed > boundaries_sorted.bed
sort -k1,1 -k2,2n intergenic.bed > intergenic_sorted.bed
sort -k1,1 -k2,2n exons.bed > exons_sorted.bed
sort -k1,1 -k2,2n fiveUTRs.bed > fiveUTRs_sorted.bed
sort -k1,1 -k2,2n threeUTRs.bed > threeUTRs_sorted.bed
sort -k1,1 -k2,2n onetofivekb.bed > onetofivekb_sorted.bed
echo "Done sorting"

echo "Merging bed files"
bedtools merge -i enhancers_sorted.bed > enhancers_merged.bed
bedtools merge -i promoters_sorted.bed > promoters_merged.bed
bedtools merge -i introns_sorted.bed > introns_merged.bed
bedtools merge -i boundaries_sorted.bed > boundaries_merged.bed
bedtools merge -i intergenic_sorted.bed > intergenic_merged.bed
bedtools merge -i exons_sorted.bed > exons_merged.bed
bedtools merge -i fiveUTRs_sorted.bed > fiveUTRs_merged.bed
bedtools merge -i threeUTRs_sorted.bed > threeUTRs_merged.bed
bedtools merge -i onetofivekb_sorted.bed > onetofivekb_merged.bed
echo "Done merging"

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
echo "Done enriching for CpG annotations"

echo "Removing temporary files"
rm *_sorted.bed
rm *_merged.bed
echo "Done removing temporary files"

echo "Done enriching for genic annotations"

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
