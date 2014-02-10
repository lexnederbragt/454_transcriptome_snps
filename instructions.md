###Principles

1 map all cDNA reads from one species to the zebrafinch transcriptome  
2 use the consensus called on the reads mapped as reference for mapping all cDNA reads from *both* species  
3 find variant bases which are only called in one of the species, i.e., all reads from one species show a different base than from the other species  
4 repeat steps 1-3, but starting with the reads from the other species

####Mapping
Used as reference genome all cDNAs from zebrafinch, *Taeniopygia_guttata*, file Taeniopygia_guttata.taeGut3.2.4.62.cdna.biomart.fa downloaded from ensemble release 62 through biomart.

For mapping, Newbler version 2.5.3 was used.

`DOM` stands for house sparrow, *P. domesticus*
`HISP` stands for Spanish sparrows, *P. hispaniolensis* 

####Prerequiste files

* MID config files, describing the tags used for multiplexing, in the format needed for Newbler

```
DOM_MIDConfig.parse
HISP_MIDConfig.parse
JoMID_MIDConfig.parse
```

* 454 SFF files, see the Hermansen et al paper for ENA/SRA accession number(s). Make one file with all sff files (the MID config files are used to select reads for the mappings)

```
all_reads.sff
```

* File mapping all reads to samples

```
readID_sample.tsv
```

####Step 1: DOM reads to zebrafinch transcripts
```
PROJECT=110621_DOM_vs_Tg_ml90%mi95
newMapping -cdna $PROJECT
addRun -mcf DOM_MIDConfig.parse $PROJECT JoMID@/data/all_reads.sff
setRef -cref $PROJECT Taeniopygia_guttata.taeGut3.2.4.62.cdna.biomart.fa
runProject -ml 90% -mi 95  $PROJECT
```

File used for next step: `454AllContigs.fna`, containing the consensus sequences of the mapped reads.

####Step 2: mapping of all DOM + HISP reads to DOM-mapped consensus contigs

**Adjust names of consensus sequences**

Consensus sequences from first mapping become references for second mapping. This step changes the names of sequences from first mapping to include name of the species, e.g.

`>contig00001  ENSTGUG00000000018, 1..2846  length=2845   numreads=89`

to

`>DOM00001_ENSTGUG00000000018_1..2846`

```
cd 110621_DOM_vs_Tg_ml90%mi95/mapping/
cat 454AllContigs.fna|awk 'BEGIN{OFS="_"}{if ($0 ~ />/) {gsub (/contig/,"DOM",$1); print $1,substr($2,1,18),$3} else {print $0}}' >DOM_vs_Tg_AllContigs.fna
```

This file is part of the repository.

**Second mapping**

```
PROJECT=110622_All_vs_DOM_mapped_ctgs_ml90%mi95
newMapping -cdna $PROJECT
addRun -mcf JoMID_MIDConfig.parse $PROJECT JoMID@/data/all_reads.sff
setRef -cref $PROJECT 110621_DOM_vs_Tg_ml90%mi95/mapping/DOM_vs_Tg_AllContigs.fna
runProject -ud -sio -ml 90% -mi 95 $PROJECT
```

####Step 3: find variant bases which are only called in one of the species

The `454HCDiffs.txt` contains the differences between the reads and the references, and lists all the mapped reads, divided between reads that map with (a) difference(s), or show the identical base(s) at the variable positions.  
To summarise read counts for each sample (HISP* and DOM*), use the  `parse454Diffs.pl` script. **NOTE**, the path to the `readID_sample.tsv` file is hardcoded in the script (sorry…)

```
cd 110830_All_vs_HISP_mapped_ctgs_ml90%mi95/mapping/
path/to/scripts/parse454Diffs.pl 454HCDiffs.txt >HCDiffs_with_sample.tsv
```

The `HCDiffs_with_sample.tsv` is part of the repository and file has the following columns:

* Reference Accno (e.g. DOM00003_ENSTGUG00000000025_22..1412, i.e. consensus contig from with DOM reads mapped to zebrafinch transcript ENSTGUG00000000025 position 22 to 1412)
* Start position of variant in the reference contig
* End position of variant in the reference contig
* Reference nucleotide(s)
* Variant nucleotide(s)
* Total read depth
* Variant frequency (frequency of reads with variant)
* 12 pairs of columns with Ref and Var read counts for all six DOM, and all six HISP read datasets. E.g. Ref DOM3 0 Var DOM3 4 would mean no reads from sample DOM3 mapped to this place with the reference base, and 4 reads from the same sample mapped with the variant base.
* four columns (Ref DOM_count, Var DOM_count, Ref HISP_count, Var HISP_count) giving the total counts for reads of the DOM samples, and HISP samples, showing the Ref or Var base, respectively

The logic is that the best candidate SNPs that distinguish between DOM and HISP are present in all or most samples of one species, but never seen in any of the other. Practically, this means those SNPs that:

* **Case 1:** show a zero for both Ref DOM_count and Var HISP_count, but not zero for Var DOM_count and Ref HISP_count,
* **Case 2:**  show a zero for both Var DOM_count and Ref HISP_count, but not zero for Ref DOM_count and Var HISP_count

A series of awk commands was used to generate a file with these SNPs.

First, copy the header and add one more column header:

```
head -2 HCDiffs_with_sample.tsv |awk 'BEGIN{OFS=FS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$32,$33,$34,$35,"sum_samples"}' >relevant_SNPs.tsv
```
Then, select SNPs:

**Case 1**

```
cat HCDiffs_with_sample.tsv|awk 'BEGIN{OFS=FS="\t"}$33==0&&$34==0 {print $1, $2,$3,$4,$5,$6,$7,$32,$33,$34,$35, $32+$33+$34+$35}' >>relevant_SNPs.tsv
```

**Case 2**

```
cat HCDiffs_with_sample.tsv|awk 'BEGIN{OFS=FS="\t"}$32==0&&$35==0 {print $1, $2,$3,$4,$5,$6,$7,$32,$33,$34,$35, $32+$33+$34+$35}' >>relevant_SNPs.tsv
```

The more samples have reads mapping to the SNP region, the better the candidate. So, we sort by number of samples having reads mapped to the SNP region:

```
mv relevant_SNPs.tsv temp.tsv
sort -nr -k12,12 temp.tsv >relevant_SNPs.tsv
rm temp.tsv
```

Then manually put the headers on top again.

Splitting the first cell, e.g. from

`>HISP03314_ENSTGUG00000004885_34..2952`

to

`contig03314     ENSTGUG00000004885 34 2952`

```
cat relevant_SNPs.tsv | perl -e 'while (defined($line=<>)){$line =~ />HISP(\d+)_(.*)_(\d+)..(\d+)(.*)/;print join "\t", ("contig".$1,$2,$3,$4,$5); print "\n"}' >relevant_SNPs_2.tsv
```

Manually:
* adjust the headers.
* remove empty fifth column
* change column order (contig number to fourth column)

**Adding flanking sequence **  
The `get_SNP_flanks.pl` perl script adds the flanking sequences (max 100 bp on either side) and the description of the zebrafinch gene the reads mapped to.
**NOTE** the path to the Taeniopygia_guttata.taeGut3.2.4.62.cdna.biomart.fa files is hardcoded in the script (sorry…).

```
path/to/scripts/get_SNP_flanks.pl 110622_All_vs_DOM_mapped_ctgs_ml90%mi95 relevant_SNPs_2.tsv >relevant_SNPs_final.tsv
```

**Column orders in the `relevant_SNPs_final.tsv` file**

Column name |Description
------------|--------
DOM/HISP_mapping_contig_number|from 110829_HISP_vs_Tg_ml90%mi95
DOM/HISP_contig_start|start position of the variant in the contig
DOM/HISP_contig_end|end position
left_flank|max 100 bp to the left of the SNP location
SNP|e.g. [T/C]
right_flank|max 100 bp to the right of the SNP location
Total_Depth|Total read depth
Var_Freq|% of reads with variant
Ref_DOM_count|Number of DOM samples showing the Ref base
Var_DOM_count|Number of DOM samples showing the Var base
Ref_HISP_count|Number of HISP samples showing the Ref base
Var_HISP_count|Number of HISP samples showing the Var base
sum_samples|Total number of either DOM or HISP samples involved (max 12)
ENSTGUG_reference|from 110829_HISP_vs_Tg_ml90%mi95
ENSTGUG_map_start|where contig starts relative to ENSTGUG
ENSTGUG_map_end|end
description|from the ENSTGUG transcript the reads mapped to (extracted from `Taeniopygia_guttata.taeGut3.2.4.62.cdna.biomart.fa` file)


The `relevant_SNPs_final.tsv` file is part of this repository.

####Step 4 repeat steps 1-3, but starting with the reads from the other species

For **Step 1**, change the `addRun` to 

```
addRun -mcf HSIP_MIDConfig.parse $PROJECT JoMID@/data/all_reads.sff
```

For **Step 2 and 3**, use `HISP` when changing the contigs names in the `454AllContigs.fna`.
