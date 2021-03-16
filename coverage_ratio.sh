input_file=$1
reference=$2
STAR --genomeDir $reference --readFilesIn $input_file --outSAMtype BAM SortedByCoordinate
echo "depth"
samtools depth Aligned.sortedByCoord.out.bam > depth
awk -v OFS="\t" '{print $1,$2-1,$2,$3}' depth > depth_inter
echo "intersect"
rm depth
bedtools intersect -wo -a depth_inter -b ~/reference/exononly > all_inter
rm depth_inter
awk -v OFS="\t" '$9=="+"{print$3,$4,$8,$9}' all_inter > allinter_pos
awk -v OFS="\t" '$9=="-"{print$3,$4,$8,$9}' all_inter > allinter_neg
rm all_inter
sort -k 3,3 -k 1,1nr allinter_pos > all_possorted
awk -v OFS="\t" '{print NR,$2,$3,$4}' all_possorted > pos_nr
sort -k 3,3 -k 1,1n allinter_neg >all_negsorted
awk -v OFS="\t" '{print NR,$2,$3,$4}' all_negsorted > neg_nr
cat pos_nr neg_nr > sorted_all
rm all_possorted
rm all_negsorted
rm pos_nr
rm neg_nr
awk '$2>0{print $0}' sorted_all > sorted_alln0
awk 'BEGIN{a="";x=1}{if (NR==1){a=$3; printf a"\t"$2} else {if ($3==a){printf"\t"$2; x++} else {a=$3;print"";printf a"\t"$2}}}' sorted_alln0 > mapping_mat
awk -v OFS="\t" '{print $1,NF}' mapping_mat > for_join
#resulting for_join files of each cell can be combined to generate a matrix including mapped bases of all transcriptsd.
