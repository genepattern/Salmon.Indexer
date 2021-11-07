#!/bin/bash
# Using getopt
# Adaped from https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh
set -e

trap abort ERR PROF
abort()
{
rm -rf $outfolder

    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

threads=1
kmer=31

echo "****************"
echo "*** getDecoy ***"
echo "****************"
while getopts ":a:t:g:k:m:f:l:j:" opt; do
    case $opt in
        a)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            gtffile=`realpath $OPTARG`
            echo "-a <Annotation GTF file> = $gtffile"
          elif [[ $OPTARG =~ ^-. ]];then
            gtffile=""
            let OPTIND=$OPTIND-1
          else
            gtffile=`realpath $OPTARG`
            echo "-t <-a <Annotation GTF file> = $gtffile"
          fi
            ;;
        t)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            txpfile=`realpath $OPTARG`
            echo "-t <Transcriptome fasta> = $txpfile"
          elif [[ $OPTARG =~ ^-. ]];then
            txpfile=""
            let OPTIND=$OPTIND-1
          else
            txpfile=`realpath $OPTARG`
            echo "-t <Transcriptome fasta> = $txpfile"
          fi
            ;;
        g)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            genomefile=`realpath $OPTARG`
            echo "-g <Genome fasta = $genomefile"
          elif [[ $OPTARG =~ ^-. ]];then
            genomefile=""
            let OPTIND=$OPTIND-1
          else
            genomefile=`realpath $OPTARG`
            echo "-t <Genome fasta> = $genomefile"
          fi
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer size> = $kmer"
            ;;
        m)
            mode="$OPTARG"
            echo "Building $mode decoy index"
            ;;
        f)
            format="$OPTARG"
            ;;
        l)
            libdir="$OPTARG"
            echo "-l = $libdir"
            bedtools=`realpath ${libdir}bedtools.static.binary-2.30.0`
            mashmap=`realpath ${libdir}mashmap-2.0`
            gffread=`realpath ${libdir}gffread-0.12.7`
            ;;
        j)
            threads="$OPTARG"
            echo "-j <Concurrency level> = $threads"
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            abort
            ;;
    esac
done

if [ ! -z "$gtffile" ]; then
 outfolder=$(basename $gtffile)
 outfolder=${outfolder/%.gz}
 outfolder=${outfolder/%.GTF}
 outfolder=${outfolder/%.gtf}
elif [ -z "$gtffile" ] && [ ! -z "$txpfile" ]; then
 outfolder=$(basename $txpfile)
 outfolder=${outfolder/%.gz}
 outfolder=${outfolder/%.fasta}
 outfolder=${outfolder/%.fa}
 indexable_fasta=$txpfile
else
echo "Either a gtf file and a genome fasta, or a transcriptome fasta must be provided" 
abort
fi

## Construct Transcripts fasta from GTF and Genome Fasta
if [ -z "$txpfile" ] && [ ! -z "$genomefile" ];	then

echo "Transcriptome fasta not provided, constructing one from GTF and Genome fasta..."
echo "This fasta file will be included in the Job Results."

 txpfile=$(basename $gtffile)
 txpfile=${txpfile/%.GTF.gz}
 txpfile=${txpfile/%.gtf.gz}
 txpfile=$txpfile.fa
gunzip -ck $genomefile > genome
gunzip -ck $gtffile > transcriptome

#gffread version

$gffread transcriptome -w $txpfile -g genome -G -L --w-nocds

#$bedtools version
#awk '$3 ~ /^(transcript)$/' transcriptome > tx_transcriptome
#rm -rf transcriptome
#$bedtools getfasta -name -split -fi genome -bed tx_transcriptome -fo $txpfile #Transcripts aren't annoated with transcript name

gzip $txpfile
txpfile=$txpfile.gz
rm -rf genome genome.fai transcriptome #tx_transcriptome

elif [ -z "$txpfile" ] && [ -z "$genomefile" ]; then
echo "Unable to construct a transcriptome fasta without a genome fasta provided."
abort
fi

## Construct full decoy index
if [ $mode == "full" ] && [ ! -z "$genomefile" ]; then
echo "Building full decoy index (using entire genome of the organism as the decoy sequence)..."
# extracting the names of the decoys
echo "Extracting decoy sequence ids"
grep "^>" <(gunzip -c $genomefile) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
rm decoys.txt.bak

# concatenating genome decoy to transcriptome
echo "Making gentrome"
cat $txpfile $genomefile > gentrome.fa.gz

indexable_fasta=gentrome.fa.gz

elif [ $mode == "partial" ] && [ ! -z "$genomefile" ]; then
echo "Building partial decoy index (transcriptome used to hard-mask the genome for decoy sequence)..."

# extracting all the exonic features to mask
echo "Extracting exonic features from the gtf"
awk -v OFS='\t' '{if ($3=="exon") {print $1,$4,$5}}' <(gunzip -ck $gtffile) > exons.bed

# masking the exonic regions from the genome
echo "Masking the genome fasta"
$bedtools maskfasta -fi <(gunzip -ck $genomefile) -bed exons.bed -fo reference.masked.genome.fa

# aligning the transcriptome to the masked genome
echo "Aligning transcriptome to genome"
$mashmap -r reference.masked.genome.fa -q $txpfile -t $threads --pi 80 -s 500

# extracting the bed files from the reported alignment
echo "Extracting intervals from mashmap alignments"
awk -v OFS='\t' '{print $6,$8,$9}' mashmap.out | sort -k1,1 -k2,2n - > genome_found.sorted.bed

# merging the reported intervals
echo "Merging the intervals"
$bedtools merge -i genome_found.sorted.bed > genome_found_merged.bed

# extracting relevant sequence from the genome
echo "Extracting sequences from the genome"
$bedtools getfasta -fi reference.masked.genome.fa -bed genome_found_merged.bed -fo genome_found.fa

# concatenating the sequence at per chromsome level to extract decoy sequences
echo "Concatenating to get decoy sequences"
awk '{a=$0; getline;split(a, b, ":");  r[b[1]] = r[b[1]]""$0} END { for (k in r) { print k"\n"r[k] } }' genome_found.fa > decoy.fa

# concatenating decoys to transcriptome
echo "Making gentrome"
cat $txpfile decoy.fa > gentrome.fa
gzip gentrome.fa.gz

# extracting the names of the decoys
echo "Extracting decoy sequence ids"
grep ">" decoy.fa | awk '{print substr($1,2); }' > decoys.txt

# removing extra files
echo "Removing temporary files"
rm exons.bed reference.masked.genome.fa mashmap.out genome_found.sorted.bed genome_found_merged.bed genome_found.fa decoy.fa reference.masked.genome.fa.fai

indexable_fasta=gentrome.fa

fi

echo "Indexing gentrome on kmer size: $kmer"
mkdir -p $outfolder

params=()
[[ $format == "TRUE" ]] && params+=(--gencode)

# If a genome is available, then index with decoys
if [ ! -z "$genomefile" ]; then
salmon index \
      --no-version-check \
      -t $indexable_fasta \
      -d decoys.txt \
      "${params[@]}" \
      -p $threads \
      -i $outfolder \
      -k $kmer ;

# removing extra files
rm decoys.txt gentrome.fa.gz

echo "Compressing Index"
tar -czvf ${outfolder}.k${kmer}.salmon_${mode}_decoy_index.tar.gz -C $outfolder .

else
salmon index \
      --no-version-check \
      -t $indexable_fasta \
      "${params[@]}" \
      -p $threads \
      -i $outfolder \
      -k $kmer ;

echo "Compressing Index"
tar -czvf ${outfolder}.k${kmer}.salmon_index.tar.gz -C $outfolder .

fi


echo "Cleaning up"
rm -rf $outfolder

trap : 0
echo >&2 '
***********************
*** DONE Processing ***
***********************
'
