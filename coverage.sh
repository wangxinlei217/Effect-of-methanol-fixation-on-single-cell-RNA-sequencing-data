#####calculate depth file
#!/bin/bash
#####need to generate a exononly file to extract bases belong to exons.
#####length files are made by seprating transcripts into ten groups with equal amount.

input_file=$1
lengthfile=$2
out_dir=$3/
reference=$4

mkdir $out_dir
STAR --genomeDir $4/STAR_reference --readFilesCommand zcat --readFilesIn $input_file --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $3


echo "====calculate depth for each transcript===="
samtools depth $out_dir/Aligned.sortedByCoord.out.bam > $out_dir/depth
echo "=====bedtools intersecting====="
awk -v OFS="\t" '{print $1,$2-1,$2,$3}' $out_dir/depth > $out_dir/depth_inter
bedtools intersect -wo -a $out_dir/depth_inter -b $4/exononly > $out_dir/all_inter
echo "====seprating inter file======="
awk -v OFS="\t" '$9=="+"{print$3,$4,$8,$9}' $out_dir/all_inter > $out_dir/allinter_pos
awk -v OFS="\t" '$9=="-"{print$3,$4,$8,$9}' $out_dir/all_inter > $out_dir/allinter_neg


j=1
while [ $j -le 178210 ];do
echo "=====getting inter within one range===="


awk 'NR==FNR{a[$1];next} $3 in a{print $0}' $lengthfile/lengtha$j $out_dir/allinter_pos >$out_dir/length_pos$j
awk 'NR==FNR{a[$1];next} $3 in a{print $0}' $lengthfile/lengtha$j $out_dir/allinter_neg >$out_dir/length_neg$j

echo "====sort and make depth matrix====="
sort -k 3,3 -k 1,1nr $out_dir/length_pos$j > $out_dir/length_possorted$j
awk -v OFS="\t" '{print NR,$2,$3,$4}' $out_dir/length_possorted$j > $out_dir/pos_nr$j

sort -k 3,3 -k 1,1n $out_dir/length_neg$j > $out_dir/length_negsorted$j
awk -v OFS="\t" '{print NR,$2,$3,$4}' $out_dir/length_negsorted$j > $out_dir/neg_nr$j
cat $out_dir/pos_nr$j $out_dir/neg_nr$j > $out_dir/sorted_length$j

awk 'BEGIN{a="";x=1}{if (NR==1){a=$3; printf a"\t"$2} else {if ($3==a){printf"\t"$2; x++} else {a=$3;print"";printf a"\t"$2}}}' $out_dir/sorted_length$j > $out_dir/$i/mapping_mat$j

awk '{$1="";print $0}' $out_dir/mapping_mat$j > $out_dir/mat$j
awk '{for (i=1;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print " "i"  ""\t"sum[i];}' $out_dir/mat$j | sort -k 1,1n > $out_dir/for_join$j



j=$((j+19801))

done

rm $out_dir/length_*
rm $out_dir/allinter*
rm $out_dir/neg_*
rm $out_dir/pos_*





###combine transcript length file from all cells into one

#!/bin/bash
##resulting lengthmat files can be load into R for visualization.


output_dir=$1

j=1
while [ $j -le 178210 ];do


pwd .
find . -name "for_join$j" -exec paste {} + > $output_dir/joined$j
find . -name "for_join$j" -type f | wc -l
awk -v OFS="\t" '{for(i=1;i<=NF;i+=2) $i="" }1' $output_dir/joined$j > $output_dir/lengthmat$j

awk '{print NR"\t"$0}' $output_dir/lengthmat$j > $output_dir/lengthmat$j

j=$((j+19801))

done
