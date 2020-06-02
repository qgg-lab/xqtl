#! /bin/bash
# ============================================================

# pipeline to perform analysis of a flyland-like experiment
# 1. alignment using bwa mem
# 2. gatk postprocessing
# 3. count alleles and perform statistics
# store only the final recal.bam files and allele counts
# all logs should be kept

# require bwa, gatk, perl (some modules too), and R
# need the perl script (requires a fisher's exact test package) to process allele counts
# and the R script to do statistics

# ============================================================

# get command line options
args=`getopt -o "b:,m:,g:,p:,e:,x:,r:,c:,u:,f:,s:,q:,i:,v:,t:,n:,l:,h:,d:,y:,z:" -l "bwa:,sam:,gatk:,picard:,bed:,index:,ref:,conv:,pileup:,filter:,stat:,fastq:,indel:,var:,tmp:,ncpu:,lg:,hg:,site:,nchr1:,nchr2:" -- "$@"`
echo "command arguments given: $args"
eval set -- "$args"

# parse arguments
while true;
do
  case $1 in
    -b|--bwa)
       bwa=$2
       shift 2;;

    -m|--sam)
       sam=$2
       shift 2;;

    -g|--gatk)
       gatk=$2
       shift 2;;

    -p|--picard)
       picard=$2
       shift 2;;

    -e|--bed)
       bed=$2
       shift 2;;

    -x|--index)
       index=$2
       shift 2;;

    -r|--ref)
       ref=$2
       shift 2;;

    -c|--conv)
       conv=$2
       shift 2;;

    -u|--pileup)
       pileup=$2
       shift 2;;

    -f|--filter)
       filter=$2
       shift 2;;

    -s|--stat)
       stat=$2
       shift 2;;

    -q|--fastq)
       fastq=$2
       shift 2;;

    -i|--indel)
       indel=$2
       shift 2;;

    -v|--var)
       var=$2
       shift 2;;

    -t|--tmp)
       tmp=$2
       shift 2;;

    -n|--ncpu)
       ncpu=$2
       shift 2;;

    -l|--lg)
       lg=$2
       shift 2;;

    -h|--hg)
       hg=$2
       shift 2;;

    -d|--site)
       site=$2
       shift 2;;

    -y|--nchr1)
       nchr1=$2
       shift 2;;

    -z|--nchr2)
       nchr2=$2
       shift 2;;

    --)
       shift
       break;;
  esac
done

# check if files and programs exist
# bwa
bwacheck=`command -v $bwa`
samcheck=`command -v $sam`
bedcheck=`command -v $bed`

if [[ $bwacheck ]]
then
  echo "using $bwacheck"
else
  echo "cannot find bwa"
  exit 1
fi

if [[ $samcheck ]]
then
  echo "using $samcheck"
else
  echo "cannot find samtools"
  exit 1
fi

if [[ $bedcheck ]]
then
  echo "using $bedcheck"
else
  echo "cannot find bedtools"
  exit 1
fi

# gatk
if [[ -e $gatk ]]
then
  echo "using $gatk"
else
  echo "cannot find $gatk"
  exit 1
fi

# picard
if [[ -e $picard/MarkDuplicates.jar ]]
then
  echo "using $picard/MarkDuplicates.jar"
else
  echo "cannot find $picard/MarkDuplicates.jar"
  exit 1
fi

# bwa index
if [[ -e $index.amb && -e $index.ann && -e $index.bwt && -e $index.pac && -e $index.sa ]]
then
  echo "using $index.*"
else
  echo "cannot find one of $index.amb .ann .bwt .pac .sa"
  exit 1
fi

# reference sequence
if [[ -e $ref.fasta && -e $ref.fasta.fai && -e $ref.dict ]]
then
  echo "using $ref.*"
else
  echo "cannot find one of $ref.fasta .fasta.fai .dict"
  exit 1
fi

# site to pileup
if [[ -e $site ]]
then
    echo "using $site"
else
    echo "cannot find $site"
    exit 1
fi

# filter script
if [[ -e $filter ]]
then
  echo "using $filter"
else
  echo "cannot find $filter"
  exit 1
fi

# statistics script
if [[ -e $stat ]]
then
  echo "using $stat"
else
  echo "cannot find $stat"
  exit 1
fi

# check chromosome number
if [[ $nchr1 =~ ^[0-9]+$ ]]
then
    echo "first sample $nchr1 chromosomes"
else
    echo "nchr1 must be a number"
    exit 1
fi

if [[ $nchr2 =~ ^[0-9]+$ ]]
then
    echo "second sample $nchr2 chromosomes"
else
    echo "nchr1 must be a number"
    exit 1
fi

# check if all fastq files exist
# the fastq file needs to follow this format
# rgid sm ln pu pair(?) illumina(I)/sanger(S) file1 file2
if [[ -e $fastq ]]
then
  echo "testing file existence for fastq files"
  while read line
  do
    pair=`echo $line | awk '{print $5}'`
    file1=`echo $line | awk '{print $7}'`
#    file1path=`eval realpath $file1`
    if [[ ! -e $file1 ]]
    then
      echo "cannot find $file1"
      exit 1
    fi
    if [[ $pair = "pair" ]]
    then
      file2=`echo $line | awk '{print $8}'`
#      file2path=`eval realpath $file2`
      if [[ ! -e $file2 ]]
      then
        echo "cannot find $file2"
        exit 1
      fi
    fi
  done < $fastq
else
  echo "cannot find the fastq input file list"
  exit 1
fi

# indel vcf
if [[ -e $indel ]]
then
  echo "using $indel"
else
  echo "cannot find $indel"
  exit 1
fi

# know variants
if [[ -e $var ]]
then
  echo "using $var"
else
  echo "cannot find $var"
  exit 1
fi

# check if fisher's exact test perl module exists
perl -e 'use Text::NSP::Measures::2D::Fisher::twotailed;'

if [ $? -ne 0 ]
then
  echo "cannot find perl module Text::NSP::Measures::2D::Fisher::twotailed"
  exit 1
fi

# ============================================================
# main program

if [[ ! -e "map/" ]]
then
  echo "map does not exist, making one"
  mkdir map
else
  echo "map already exists"
  exit 1
fi

# 1. make alignment

while read line
do

  rgid=`echo $line | awk '{print $1}'`
  sm=`echo $line | awk '{print $2}'`
  lb=`echo $line | awk '{print $3}'`
  pu=`echo $line | awk '{print $4}'`
  pair=`echo $line | awk '{print $5}'`
  file1=`echo $line | awk '{print $7}'`
  file2=""
  qual=`echo $line | awk '{print $6}'`

  if [[ $qual == "I" ]]
  then
    gunzip -c $file1 | perl $conv illumina2std > map/file1.fq
  else
    gunzip -c $file1 > map/file1.fq
  fi

  if [[ $pair == "pair" ]]
  then
    file2=`echo $line | awk '{print $8}'`
    if [[ $qual == "I" ]]
    then
      gunzip -c $file2 | perl $conv illumina2std > map/file2.fq
    else
      gunzip -c $file2 > map/file2.fq
    fi
    $bwa mem -M -t $ncpu -R "@RG\tID:$rgid\tPL:illumina\tLB:$lb\tSM:$sm\tPU:$pu" $index map/file1.fq map/file2.fq 2> map/$rgid.mem.log | $sam view -bS -t $ref.fai - 2> map/$rgid.samview.log | $sam sort -m 20000000000 - map/$rgid.map 2> map/$rgid.samsort.log
    rm map/file1.fq map/file2.fq    
  else
    $bwa mem -t $ncpu -R "@RG\tID:$rgid\tPL:illumina\tLB:$lb\tSM:$sm\tPU:$pu" $index map/file1.fq 2> map/$rgid.mem.log | $sam view -bS -t $ref.fai - 2> map/$rgid.samview.log | $sam sort -m 20000000000 - map/$rgid.map 2> map/$rgid.samsort.log
    rm map/file1.fq
  fi

  # echo -e "@RG\tID:$rgid\tPL:illumina\tLB:$lb\tSM:$sm\tPU:$pu" >> map/all.sample.sam.header

done < $fastq

# 2. indel realignment

for rg in $(awk '{print $1}' $fastq)
do
  $sam index map/$rg.map.bam
  java -jar $gatk \
    -I map/$rg.map.bam \
    -R $ref.fasta \
    -T RealignerTargetCreator \
    -nt $ncpu \
    -o map/$rg.indel.intervals > map/$rg.indel.target.log 2>&1
done

# merge indel intervals
cat map/*.indel.intervals | sed 's/:/\t/' | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[1] =~ m/-/) { @pos = split /-/, $line[1]; print $line[0], "\t", $pos[0], "\t", $pos[1], "\n"; } else { print $line[0], "\t", $line[1], "\t", $line[1], "\n"; }' | sort -k1,1 -k2,2n -k3,3n | $bed merge -i stdin | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[2] > $line[1]) { print $line[0], ":", $line[1], "-", $line[2], "\n"; } else { print $line[0], ":", $line[1], "\n"; }' > map/combine.indel.intervals

for rg in $(awk '{print $1}' $fastq)
do
  java -jar $gatk \
    -I map/$rg.map.bam \
    -R $ref.fasta \
    -T IndelRealigner \
    -targetIntervals map/combine.indel.intervals \
    -known $indel \
    -o map/$rg.realn.bam > map/$rg.indel.realn.log 2>&1
done

# 3. remove duplicates

for rg in $(awk '{print $1}' $fastq)
do
  java -jar $picard/MarkDuplicates.jar \
    input=map/$rg.realn.bam \
    output= map/$rg.realn.rmdup.bam \
    validation_stringency=LENIENT \
    metrics_file= map/$rg.rmdup.metrics \
    remove_duplicates=false assume_sorted=true \
    max_file_handles_for_read_ends_map=512 \
    create_index=true tmp_dir=$tmp \
    >  map/$rg.rmdup.log 2>&1
done

# 4. Count covariates

for rg in $(awk '{print $1}' $fastq)
do
  java -Djava.io.tmpdir=$tmp -jar $gatk \
    -T BaseRecalibrator \
    -I map/$rg.realn.rmdup.bam \
    -R $ref.fasta \
    -knownSites $var \
    -nct $ncpu \
    -o map/$rg.recal.data.grp > map/$rg.countCovariate.log 2>&1
done

# 5. Recalibration

for rg in $(awk '{print $1}' $fastq)
do
  java -Djava.io.tmpdir=$tmp -jar $gatk \
    -T PrintReads \
    -R $ref.fasta \
    -I map/$rg.realn.rmdup.bam \
    -BQSR map/$rg.recal.data.grp \
    -nct $ncpu \
    -o map/$rg.recal.bam > map/$rg.recal.log 2>&1  
done

# 6. delete intermediate files
rm map/*.map.bam
rm map/*.realn.b*
rm map/*.realn.rmdup.b*

# 7. Pile up and filter

$sam mpileup -Q 1 -q 13 -d 10000 -l $site <(awk '$2 == "'$lg'"' $fastq | awk '{print "map/"$1".recal.bam"}' | xargs $sam merge -) <(awk '$2 == "'$hg'"' $fastq | awk '{print "map/"$1".recal.bam"}' | xargs $sam merge -) 2> map/pileup.log | perl $pileup -minQ 13 -track map/pileup.track.out -ref $ref.fasta -mono 0 2> map/cov.log | perl $filter -minGB 10 -minTC 10 2> map/filter.log | awk '($6 == "AF" || $6 == "PASS") && ($8 == "AF" || $8 == "PASS")' > map/pileup.out

# 8. BSA test
Rscript $stat map/pileup.out $nchr1 $nchr2 map/bsa.out > map/bsa.log 2>&1 
