# immunoseq2airr

### v 1.2.0
#### Jamie Heather, 2021, MGH

[![DOI](https://zenodo.org/badge/258680195.svg)](https://zenodo.org/badge/latestdoi/258680195)

This script takes **v2** data from Adaptive Biotech's immunoSEQ/immuneACCESS platform, and converts them into the standardised [AIRR Community rearrangement format](https://docs.airr-community.org/en/latest/) proposed in [Vander Heiden *et al.*](https://doi.org/10.3389/fimmu.2018.02206) The goal of this effort is to make data produced by [Adaptive Biotech](https://www.adaptivebiotech.com/) compatible with ongoing standardisation efforts, and increase data re-use and comparability. To that end, during the format conversion TCR gene names are also converted into their proper HUGO-approved IMGT equivalents[.](http://jamimmunology.blogspot.com/2018/09/the-problem-with-adaptive-tcr-data-nomenclature.html)

As this is being done on the final output as provided in the '**v2**' export option, and not using any intermediate format, there will be a number of missing fields (e.g. any of the alignment descriptor fields required to be present by the standard). It is also only built for TCR sequence data currently, although it should be theoretically adaptable to immunoglobulins.

It currently also only runs on human TRA and TRB (alpha and beta) chain files, although it should be readily adaptable to the other loci. Users would just need to edit the '**adaptive_v_convert**' dictionary to include mappings for the other TCR genes for which the Adaptive renaming has erroneously implied gene-family status.

Please note that I do not have any affiliation with Adaptive Biotech (beyond being an occasional customer). ImmunoSEQ® and immuneACCESS are the property of Adaptive Biotech (Seattle, US). I also do not presume to speak for the AIRR-seq community, and this script was not produced as part of their efforts - merely in line with them.

This script doesn't use any packages that are likely to be absent from a typical installation. 

### Arguments

* `-i`: path to input v2 exported Adaptive data file (the only required option)
* `-n`: name to use as prefix in the sequence_id field 
* `-o`: output file name - if not provided the input filename with an added '-airr' will be used 
* `-lz`: how many leading zeroes to use in the sequence_id field (default = 8)
* `-z`: option to toggle on gzip compression of output file
* `-or`: toggle on ignoring orphons (see below)
* `-d`: prevent ambiguous D gene calls (see below)
* `-nd`: prevent all D gene calls
* `-c`: TCR chain (Alpha/A/TRA/a or Beta/B/TRB/b) (default = b)
* `-a`: allow ambiguity in allele level calls (see below)
* `-af`: abundance filter (float), to only keep rearrangements with a certain frequency
* `-pf`: productivity filter, to only keep in-frame rearrangements
* `-sa`: strip alleles, to keep only gene-level information
* `-mf`: motif filter, to discard TCRs with CDR3 junctions not bound by C and F (or W for TRA)
* `-jlf`: junction length filter, to discard TCRs with CDR3 junctions shorter than a given value
* `-p`: path to a tab-separated parameter file, if applying to a non-standard format (see below)

### Example usage

In order to run these examples please download the following files from Adaptive's immuneACCESS database using the option 'Export > Export Sample (v2)' :

* **D149_SP_CD8_NAIVE** (TRA)
    * From the project ['Long-term maintenance of human naive T cells through in situ homeostasis in lymphoid tissue sites'](https://doi.org/10.21417/B7G019) 
    * Data from [Thorne *et al.*](https://doi.org/10.1126/sciimmunol.aah6506)
* **TRA-D002-031-PBMC-Unsorted** (TRB)
    * From the project ['Origin and evolution of the T-cell repertoire after posttransplantation cyclophosphamide'](https://doi.org/10.21417/B75P4J)
    * Data from [Kanakry *et al.*](https://doi.org/10.1172/jci.insight.86252)
    
Unzip these in to the same directory as the immunoseq2airr.py script to be able to run the following. Then we'll focus on the TRB sample first:

```bash
# A simple straightfowrward TRB run
python immunoseq2airr.py -i D149_SP_CD8_NAIVE.tsv
```

This should convert the input format that looks like this (truncated):

|                                        nucleotide                                       |     aminoAcid     | count |  frequencyCount (%)  | cdr3Length |  vMaxResolved | vFamilyName |  vGeneName | vGeneAllele | vFamilyTies |     vGeneNameTies     | … |
|:---------------------------------------------------------------------------------------:|:-----------------:|:-----------------------:|:--------------------:|:----------:|:-------------:|:-----------:|:----------:|:-----------:|:-----------:|:---------------------:|:---:|
| AGCCCTCAGAACCCAGGG... |                   |            1            | 2.00812890581072E-05 |     46     |   unresolved  |   TCRBV12   | unresolved |             |             | TCRBV12-03,TCRBV12-04 |   |
| CACCTACACACCCTGCAG... |   CASSVTGNTEAFF   |            1            | 1.40569023406751E-05 |     39     | TCRBV04-02\*01 |   TCRBV04   | TCRBV04-02 |      1      |             |                       |   |
| CTGCAGCCAGAAGACTCG... | CASSQDRELAGGWTQYF |            1            | 1.80731601522965E-05 |     51     | TCRBV04-03\*01 |   TCRBV04   | TCRBV04-03 |      1      |             |                       |   |
...

Into this:

| sequence_id                | sequence                                                                                | v_call            | d_call            | j_call     | junction_aa       | duplicate_count | rev_comp | productive | ... |
|----------------------------|-----------------------------------------------------------------------------------------|-------------------|-------------------|------------|-------------------|-----------------|----------|------------|-----|
| D149_SP_CD8_NAIVE⎮00000001 | AGCCCTCAGAACCCAGGG... | TRBV12-3,TRBV12-4 | TRBD2\*01,TRBD2\*02 | TRBJ2-1\*01 |                   | 1               | F        | F          |     |
| D149_SP_CD8_NAIVE⎮00000002 | CACCTACACACCCTGCAG... | TRBV4-2\*01        | TRBD1\*01          | TRBJ1-1\*01 | CASSVTGNTEAFF     | 1               | F        | T          |     |
| D149_SP_CD8_NAIVE⎮00000003 | CTGCAGCCAGAAGACTCG... | TRBV4-3\*01        | TRBD2\*02          | TRBJ2-3\*01 | CASSQDRELAGGWTQYF | 1               | F        | T          |     |

Note that by under these default settings there may be rearrangements for which there is no given allele level information, such as this top row (listed just as 'TRBV12-3,TRBV12-4'). This is because immunoseq2air does not add any information that was not present in the input TSV. If they wish users could go through the data and reassign these calls to be listed as all possible alleles given the relevant sequence, but that exceeds the scope of this tool which aims primarily just to convert.

There are some other minor input fields one may specify:

```bash
# Change prefix in rearrangement sequence_id field
python immunoseq2airr.py -i D149_SP_CD8_NAIVE.tsv.gz -n example-prefix

# Change output file name
python immunoseq2airr.py -i D149_SP_CD8_NAIVE.tsv.gz -o example-outname

# Change output file name, compress output, increase number of leading zeroes in 
sequence_id
python immunoseq2airr.py -i D149_SP_CD8_NAIVE.tsv.gz -o example-compressed -z -lz 10
```

Users can see all of the options available by using the help flag: ```python immunoseq2airr.py -h```.

We can also choose to ignore orphon genes (i.e. those located outside the classical TCR loci) in ambiguous gene calls. This option effectively offers users to use the assumption that in the case of ambiguous gene calls between an orphon and it's in-locus paralogue the in-locus option is the most likely gen involved in a rearrangement. So 'TCRBV20-01,TCRBV20-or09_02' instead becomes just 'TRBV20-1'. Note that unambiguous orphons calls (i.e being the only gene called) are retained.

```bash
# Ignoring orphon genes
python immunoseq2airr.py -i D149_SP_CD8_NAIVE.tsv.gz -or
```

Looking at the alpha file allows us to explore some of the additional functionality.

```bash
# Repeating the simple run...
python immunoseq2airr.py -i TRA-D002-031-PBMC-Unsorted.tsv 

# This raises the following error message:
IOError: Unknown format on line 1! Cannot continue. 
        Ambiguity for D gene calls lacking allele info that is
         not resolved in either 'Gene' or 'Allele Ties' fields.
        Try re-running the script using the '-a' flag (to allow ambiguity),
         and check that the format of the output document is correct.
```

This error occurs as these files have one of the slightly-non-standard formats where unknown alleles are not specified (see 'Cautions' section). Seeing as these are alpha chain rearrangements which don't have D genes (or shouldn't!) we can just ignore all Ds with the ```-d``` flag, as it may be restricted to these irrelevant genes (which can have values entered, even when not necessarily using a non-TRD-shared TRAV).

Note that it can also be useful to override normal D gene handling, to ignore ambiguous gene calls, which effectively suppresses 'TRBD1,TRBD2' being called (which is often the case, and not very informative given that those are the two options and are frequently indistinguishable).

```bash
# Repeating the attempt but ignoring D genes
python immunoseq2airr.py -i TRA-D002-031-PBMC-Unsorted.tsv 

# Raises the following error message:
IOError: Unknown format on line 4052! Cannot continue. 
        Ambiguity for J gene calls lacking allele info that is
         not resolved in either 'Gene' or 'Allele Ties' fields.
        Try re-running the script using the '-a' flag (to allow ambiguity),
         and check that the format of the output document is correct.
```

Again we can see that the same error message is coming up, but this time for a J gene instead of a D, and for a rearragement much further into the file. As the error suggests, we can re-run the script using the ```-a``` flag, to allow these ambiguities instead of stopping, like so:

```bash
python -i immunoseq2airr.py -i TRA-D002-031-PBMC-Unsorted.tsv -nd -a
```

As we must when using the ```-a``` flag (see Cautions below), checking the output shows us a table formatted exactly like we expect. The only difference is that some fields on some rows (like the j_call for the 4,052nd data row) lack allele information (and failed to represent this in a consistent manner). 

In my experience using the ```-a``` flag as default on newly downloaded experiments can be dangerous, as occasionally missing allele level information has been a symptom of another cryptically organised file format.

#### Custom parameter files

While this script aims to cover the largest breadth of files in v2 "format", there are still occasionally some that don't fit. For circumstances where this is the case, users can work out the appropriate column conversions from the datafile and supply this information to the script in the form of a tab-separated file via the `-p` flag. Use this option to give the path of a file containing rows with a tab separated entry relating to each of the necessary fields required to be found in an Adaptive file. For convenience, an example file has been included in this repo (`emerson-parameters.tsv`) which gives an example that can be used to convert the 'HIP' prefixed files in the large [Emerson et al. Adaptive dataset](http://dx.doi.org/10.21417/B7001Z), e.g.:


```bash
# Repeating the simple run...
python immunoseq2airr.py -i HIP02805.tsv.gz -p emerson-parameters.tsv -a
```

Note that in order to cope with the Emerson data (which is quite far removed from the v2 standard) you have to use the `-a` option. ls


The fields present in the parameter file are:

| | |
:-----:|:-----:
sequence\_index|0
cdr3\_index|1
abundance\_index|4
productivity|2
vMaxResolved|10
dMaxResolved|13
jMaxResolved|16
vGeneNameTies|30
dGeneNameTies|33
jGeneNameTies|36
vGeneAlleleTies|31
dGeneAlleleTies|34
jGeneAlleleTies|37

If you generate your own be sure that the file contains only the rows present in the example file, in the same order.

### Cautions

Be warned that this script may throw unexpected errors, as it seems that there are inconsistencies between files from different experiments even if they are all exported as v2. You need to be particularly wary of using the ```-a``` 'allow_ambiguity' option, which effectively just takes whatever is in a given gene call cell even if it can't find an unambiguous entry. I suggest running the script twice in situations like this, once with and once without the ```-a``` flag; if it breaks on a very early line number, chances are that it's an uncatered-to file type. If it's a later line_count, then perhaps it's OK to proceed with allowing ambiguity (but still give your output data a once over to check!).

I've incorporated measures to cope with all of the file format heterogeneity that I've encountered, but without knowing how many possible versions of 'v2' files there are it's impossible for me to make the script anticipate all of them. A couple of examples of this behaviour include:

* Sometimes ambiguous '*MaxResolved*' calls are given as an empty field, sometimes as 'unresolved'. 
* In some file versions, ambiguity in the *MaxResolved* call is detailed in the '*GeneNameTies*' field, other times in the '*AlleleNameTies*' field (likely dependent on whether or not the gene is a member of a gene family or not, but it's just treated differently between different experiments).
* Sometimes (as in the TRA example used above) a file will just omit lists of potential alleles from both the *MaxResolved* and *GeneAlleleTies* options. In that example above, the 'TRAJ15' found on the 4052nd data row might theoretically be better (or at least more consistently) represented as 'TRAJ15\*01,TRAJ15\*02' - the only two known possibilities - however adding this information would be exceeding the remit of this script.

These inconsistencies seem to be less frequent in more recently uploaded datasets, but YMMV. 

Also note that Adaptive datasets are often liberal with what they determine to be an in-frame (i.e. potentially productive) rearrangement. This script requires there both to be an 'In' in the 'sequenceStatus' column and for there to be something present in the CDR3 field. However I've seen cases where there's an 'In' and the CDR3 field contains something that clearly just a partial CDR3 at best (e.g. 'YQLIW') which is extremely unlikely to represent a genuine fully-sequenced productive TCR, so post-processing sanity checks are advised.
