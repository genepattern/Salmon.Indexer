#!/bin/bash
# Using getopt
# Adaped from https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh
set -e

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
while getopts ":a:t:g:k:m:f:o:l:j:" opt; do
    case $opt in
        a)
            gtffile=`realpath $OPTARG`
            echo "-a <Annotation GTF file> = $gtffile"
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
            genomefile=`realpath $OPTARG`
            echo "-g <Genome fasta> = $genomefile"
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
        o)
            outfolder="$OPTARG"
            outfolder=${outfolder/%.gtf}
            echo "-o <Output files Path> = $outfolder"
            ;;
        l)
            libdir="$OPTARG"
            echo "-l = $libdir"
            bedtools=`realpath ${libdir}bedtools.static.binary-2.29.2`
            mashmap=`realpath ${libdir}mashmap-2.0`
            gffread=`realpath ${libdir}gffread-0.12.4`
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

## Construct Transcripts fasta from GTF and Genome Fasta
if [ -z "$txpfile" ]
then

echo "Transcriptome fasta not provided, constructing one from GTF and Genome fasta..."

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
fi

## Construct full decoy index
if [[ $mode == "full" ]]
	then
echo "Building full decoy index (using entire genome of the organism as the decoy sequence)..."
# extracting the names of the decoys
echo "Extracting decoy sequence ids"
grep "^>" <(gunzip -c $genomefile) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
rm decoys.txt.bak

# concatenating genome decoy to transcriptome
echo "Making gentrome"
cat $txpfile $genomefile > gentrome.fa.gz

elif [[ $mode == "partial" ]]
	then
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
gzip gentrome.fa

# extracting the names of the decoys
echo "Extracting decoy sequence ids"
grep ">" decoy.fa | awk '{print substr($1,2); }' > decoys.txt

# removing extra files
echo "Removing temporary files"
rm exons.bed reference.masked.genome.fa mashmap.out genome_found.sorted.bed genome_found_merged.bed genome_found.fa decoy.fa reference.masked.genome.fa.fai

fi

echo "Indexing gentrome on kmer size: $kmer"
mkdir -p $outfolder

params=()
[[ $format == "TRUE" ]] && params+=(--gencode)

salmon index \
      --no-version-check \
      -t gentrome.fa.gz \
      -d decoys.txt \
      "${params[@]}" \
      -p $threads \
      -i $outfolder \
      -k $kmer ;

# removing extra files
rm decoys.txt gentrome.fa.gz

echo "Compressing Index"
tar -czvf ${outfolder}.k${kmer}.salmon_${mode}_decoy_index.tar.gz -C $outfolder .

echo "Cleaning up"
rm -rf $outfolder

trap : 0
echo >&2 '
***********************
*** DONE Processing ***
***********************
'
