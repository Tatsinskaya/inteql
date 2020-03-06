#!/bin/bash
#INPUTS
inputdata='../data/original-data'
inputGTEX=$inputdata"/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
eQTLFolder=$inputdata"/GTEx_Analysis_v7_eQTL/"
TargetFinder=$inputdata"/TargetFinder/"
databases=$inputdata"/databases/"
HiCFolder=$inputdata"/GM12878_combined/5kb_resolution_intrachromosomal"
dbSNPFolder=$inputdata"/dbSNP/chr-header"
scriptfolder='/nfs/research1/zerbino/jhidalgo/inteql/scripts/'
regbuildgff=$inputdata'/homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20180925.gff'

#OUTPUTS
outputfolder="/nfs/research1/zerbino/jhidalgo/inteql/data/output/" ### MUST CONTAIN FINAL /
output01=$outputfolder"output01"
output02=$outputfolder"output02.tab"
#output03="NULL"
output04=$outputfolder"output04.csv.gz"
output05=$outputfolder"output05.csv.gz"
#output06="NULL"
output07=$outputfolder"output07.csv.gz"
output08=$outputfolder"output08.csv.gz"
output09=$outputfolder"output09_ALL.csv.gz"
perchroutput09=$outputfolder"finemap/"                             ### MUST CONTAIN FINAL /
output10=$outputfolder"output10.csv"
perchroutput10=$outputfolder"modeling/"                            ### MUST CONTAIN FINAL /
outputstdin=$outputfolder"stdin/"                                  ### MUST CONTAIN FINAL /
outputstderr=$outputfolder"stderr/"                                ### MUST CONTAIN FINAL /
backupfolder="/nfs/research1/zerbino/jhidalgo/inteql/data/backup/" ### MUST CONTAIN FINAL /

topgenes="5000"
chromosomes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

echo "Executing pre-HiC scripts"

if [ -d $outputfolder ]; then
  mkdir -p $backupfolder
  mv $outputfolder $backupfolder"output_$(date +"%y%m%d")"
  echo "Previous run files found, moved to backup/output_$(date +"%y%m%d")"
else
  echo "Previous run files not found, creating folders..."
fi
mkdir -p $perchroutput09
mkdir -p $perchroutput10
mkdir -p $outputstdin
mkdir -p $outputstderr

##SCRIPT 01 getTopGenesInd###
echo -e "############\nExecuting script 01_getTopGenesInd.py, time: $(date +"%H:%M:%S")"
bsub -J 'x01x' -M 1000 -o "$outputstdin"output01_my-stdin.txt -e "$outputstderr"output01_my-stderr.txt "python 01_getTopGenesInd.py $inputGTEX $output01 $topgenes"
bwait -w 'ended(x01x)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"
##SCRIPT 02 getTopGenesData###
echo -e "############\nExecuting script 02_getTopGenesData.py, time: $(date +"%H:%M:%S")"
bsub -J 'x02x' -M 1000 -o "$outputstdin"output02_my-stdin.txt -e "$outputstderr"output02_my-stderr.txt "python "$scriptfolder"02_getTopGenesData.py "$output01".npy $inputGTEX $output02"
bwait -w 'ended(x02x)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

##SCRIPT 03 ###
echo -e "############\nSkipping script 3, time: $(date +"%H:%M:%S")"

##SCRIPT 04 ###
echo -e "############\nExecuting script 04_getPairsIdTop.py, time: $(date +"%H:%M:%S")"
bsub -J 'x04x' -M 1000 -o "$outputstdin"output04_my-stdin.txt -e "$outputstderr"output04_my-stderr.txt "python "$scriptfolder"04_getPairsIdTop.py $eQTLFolder $output02 $output04"
bwait -w 'ended(x04x)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

##SCRIPT 05 ###
echo -e "############\nExecuting script 05_getPairsSlopeTop.py, time: $(date +"%H:%M:%S")"
bsub -J 'x05x' -M 4000 -o "$outputstdin"output05_my-stdin.txt -e "$outputstderr"output05_my-stderr.txt "python "$scriptfolder"05_getPairsSlopeTop.py $eQTLFolder $output04 $output05"
bwait -w 'ended(x05x)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

##SCRIPT 06 ###
echo -e "############\nSkipping script 6, time: $(date +"%H:%M:%S")"

##SCRIPT 07 ###
echo -e "############\nExecuting script 07_addLinear.py, time: $(date +"%H:%M:%S")"
bsub -J 'x07x' -M 6000 -o "$outputstdin"output07_my-stdin.txt -e "$outputstderr"output07_my-stderr.txt "python "$scriptfolder"07_addLinear.py $output02 $output05 $TargetFinder $output07 $regbuildgff"
bwait -w 'ended(x07x)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

##SCRIPT 08 ###
echo -e "############\nExecuting script 08_addHiC.py, time: $(date +"%H:%M:%S")"
bsub -J 'x08x' -M 60000 -o "$outputstdin"output08_my-stdin.txt -e "$outputstderr"output08_my-stderr.txt "python2.7 "$scriptfolder"08_addHiC.py $HiCFolder $output07 $output08"
bwait -w 'ended(x08x)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

###SCRIPT 09###
echo -e "############\nExecuting script 09_prepareFineMap.py, time: $(date +"%H:%M:%S")"
for i in $chromosomes; do
  echo -e "\nStarting with CHR "$i
  mkdir -p "$perchroutput09"chr_$i
  bsub -J "x09_"$i"x" -M 3000 -o $outputstdin"output09_my-stdin_chr"$i".txt" -e $outputstderr"output09_my-stderr_chr"$i".txt" python2.7 "$scriptfolder"09_prepareFineMap.py $databases $output08 $i $dbSNPFolder $perchroutput09"output09_"$i".csv.gz" $perchroutput09 $eQTLFolder"Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz"
done

for i in $chromosomes; do
  bwait -w "ended(x09_"$i"x)"
done

zcat "$perchroutput09"output09_*.csv.gz | sort -h | uniq | gzip -cf >$output09

echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

###SCRIPT 10###
for i in $chromosomes; do
  if [ -e "$perchroutput09"output09_$i.csv.gz ]; then
    echo -e "############\nExecuting script 10_modeling.py with chr $i, time: $(date +"%H:%M:%S")"
    bsub -J "x10_"$i"x" -o "$outputstdin"output10_my-stdin_chr$i.txt -e "$outputstderr"output10_my-stderr_chr$i.txt -M 3000 python2.7 "$scriptfolder"10_modeling.py "$perchroutput09"output09_$i.csv.gz $output05 $i "$perchroutput10"output10_$i.csv
  else
    echo -e "############\nSkipping script 10_modeling.py with chr $i due to missing file, time: $(date +"%H:%M:%S")"
    bsub -J "x10_"$i"x" -o /dev/null -e /dev/null "sleep 1"
  fi
done

echo -e "############\nExecuting script 10_modeling.py with All chromosomes, time: $(date +"%H:%M:%S")"
bsub -J 'x10_ALLx' -o "$outputstdin"output10_my-stdin_chrALL.txt -e "$outputstderr"output10_my-stderr_chrALL.txt -M 3000 python2.7 "$scriptfolder"10_modeling.py $output09 $output05 ALL "$perchroutput10"output10_ALL.csv

for i in $chromosomes; do
  bwait -w "ended(x10_"$i"x)"
done
bwait -w 'ended(x10_ALLx)'
echo -e "############\nFinished at time: $(date +"%H:%M:%S")"

# Provide Script 10 overall output
for i in $chromosomes; do
  if [ -e "$perchroutput09"output09_$i.csv.gz ]; then
    awk -F $'\t' -v CHR="$i" '{if(NR==2&&CHR==1) print $1 FS $2 FS $3 FS $4 FS $5 FS "Chromosome"; else if(NR>=3&&NR<=11) print $1 FS $2 FS $3 FS $4 FS $5 FS CHR}' "$perchroutput10"output10_$i.csv
  fi
done >$output10
awk -F $'\t' '{if(NR>=3&&NR<=11) print $1 FS $2 FS $3 FS $4 FS $5 FS "All"}' "$perchroutput10"output10_ALL.csv >>$output10

echo -e "############ RUN FINISHED ############"
echo -e "############ $(date) ############"

####### COUNT AND POSITION FILES####
declare -A scriptname
scriptname=([output04]="04_PairsID" [output05]="05_PairsSlope" [output07]="07_LinearData" [output08]="08_HiCData" [output09]="09_Finemap")
for a in output04 output05 output07 output08 output09; do
  echo -ne "Chromosome\t${scriptname[$a]}\n" >$outputfolder$a"_count.txt"
  for i in $chromosomes; do
    echo -ne "$i\t"
    zcat ${!a} | grep "^"$i"_" -c
  done >>$outputfolder$a"_count.txt"
  echo -ne "Total\t" >>$outputfolder$a"_count.txt"
  zcat ${!a} | wc -l >>$outputfolder$a"_count.txt"
  echo -ne "Genes\t" >>$outputfolder$a"_count.txt"
  zcat ${!a} | cut -f 2 -d "," | sort -u | wc -l >>$outputfolder$a"_count.txt"
  echo -ne "Chromosome\tPosition\tSource\n" >$outputfolder$a"_variant_pos.txt"
  zcat ${!a} | awk -F $'_' -v SOURCE=${scriptname[$a]} '{if(NR>1) print $1 "\t" $2 "\t" SOURCE}' | sort -h >>$outputfolder$a"_variant_pos.txt"
done

paste "$outputfolder"output*_count.txt | cut -f 1,2,4,6,8,10 >"$outputfolder"outputall_count.txt
cat $outputfolder"output"*"_variant_pos.txt" >$outputfolder"outputall_variantpos.txt"
