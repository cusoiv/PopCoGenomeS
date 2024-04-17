configfile=./phybreak_config.sh
source ${configfile}

#copy genomes from population to project_dir/genomes

#prepare list of genome names as pop_names
pop_infile_name=${OUTFOLDER_2}.txt
cp ${pop_infile_source}/${pop_infile_name} ${project_dir}/
cut -d$'\t' -f1 ${project_dir}/${pop_infile_name} > pop_names

#copy genomes to input folder $project_dir/genomes
while read line;do
    echo $line
    cp ${genome_source}/${line}${genome_ext} ${project_dir}/genomes
done<pop_names

#if genomes are zipped, unzip them and remove the zipped files, also create new environment variable genome_ext omitting the gz

#remove .gz from genome_ext
if [[ "$genome_ext" == *.gz ]]; then
    echo $genome_ext
    gunzip ${project_dir}/genomes/*.gz
    genome_ext="${genome_ext%.gz}"
fi
echo $genome_ext

#find ref_isolate
seqkit stats -b -N 50 ${project_dir}/genomes/*${genome_ext} > genome_info
grep -f ${project_dir}/align_and_construct_trees/pop_names genome_info > genome_info_select
sed -i 's/'"${genome_ext}"'//g' genome_info_select
ref_iso=$(awk -v max=0 '{if($NF>max){want=$1; max=$NF}}END{print want}' genome_info_select)
ref_contig=${ref_iso}_1
echo $ref_iso

#stitch ref_genome together
#head ${project_dir}/genomes/${ref_iso}${genome_ext}
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${project_dir}/genomes/${ref_iso}${genome_ext} > ${project_dir}/genomes/${ref_iso}_linear.fna
awk -F'\t' '{print $2}' ${project_dir}/genomes/${ref_iso}_linear.fna > ${project_dir}/genomes/newOH.fna
awk -v ref_contig="$ref_contig" 'BEGIN { ORS="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; print ">" ref_contig "\n" } { print }' < ${project_dir}/genomes/newOH.fna > ${project_dir}/genomes/${ref_iso}_1.fna
rm ${project_dir}/genomes/${ref_iso}.fna
rm ${project_dir}/genomes/${ref_iso}_linear.fna
rm ${project_dir}/genomes/newOH.fna
mv ${project_dir}/genomes/${ref_iso}_1.fna ${project_dir}/genomes/${ref_iso}${genome_ext}

#edit phybreak_parameter file
sed "s|.project_dir|$project_dir|g" phybreak_parameters_template.txt > phybreak_parameters.txt
sed -i "s|.pop_infile_name|$pop_infile_name|g" phybreak_parameters.txt
sed -i "s|.basename|$basename|g" phybreak_parameters.txt
sed -i "s|.ref_iso|$ref_iso|g" phybreak_parameters.txt
sed -i "s|.ref_contig|$ref_contig|g" phybreak_parameters.txt
sed -i "s|.input_contig_extension|$genome_ext|g" phybreak_parameters.txt
