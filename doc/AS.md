
## Allele specific analysis

If the ALLELE_SPECIFIC_SNP option is defined in the configuration file, HiC-Pro will run the allele specific mode.  
The ALLELE_SPECIFIC_SNP option must contain the path to a VCF file with SNPs information. The [extract_snps.py](UTILS.md) utility can be used to generate such file.  
This utility was first design to extract relevant information from VCF file provided by the [Mouse Sanger database](http://www.sanger.ac.uk/resources/mouse/genomes/). It aims at generating a VCF file of the F1 individual based on its parental genotype. For instance using *extract_snps.py* with *-r CASTEiJ* and *-a 129S1* will generate a VCF file with all F1 heterogyzote SNPs which can be used to distinguish *CASTEiJ* and *129S1* alleles.  
The idea is therefore to have a VCF will maternal/paternal haplotypes encoded as the reference/alternative SNPs information.
Phasing data, such as the ones available from the [Illumina Platinum Project](http://www.illumina.com/platinumgenomes/) can be simply used as is.


### Allele specific mapping

In allele-specific mode, the sequencing reads are first aligned on a masked reference genome for which all polymorphic sites were first N-masked.
In the current version, **this is the user responsability to generate this genome and to provide the bowtie2 indexes to HiC-Pro.** Example of how to generate such reference genome is discribed in the [FAQ](FAQ.md) section.
This masking strategy avoid systematic bias toward the reference allele, compared to standard mapping where reads with the reference allele are more likely to be mapped than the reads with non-reference alleles.

### Assignment to parental genome

Once aligned HiC-Pro browses all reads spanning a polymorphic site, locates the nucleotide at the appropriate position, and assigns the read either to the maternal or paternal allele. Reads with conflicting allele assignment or unexpected allele at polymorphic sites are discarded.  
This step generates a bam file with all reads flagged according to its parental assignment (AS flag).

* **XA:i:0** - unassigned (UA)
* **XA:i:1** - assigned to reference genome (G1)
* **XA:i:2** - assigned to alternative genome (G2)
* **XA:i:3** - conflicting (C)

### Allele specific maps

Finally, read pairs with at least one allele specific mate are used to construct allele specific contact maps.
In allele-specific mode, the *.validPairs* files will contain an additional columns with the allele status of each reads pair. Therefore, valid pairs assigned to **'G1-G1', 'G1-UA' or 'UA-G1'** are used to build the **G1** genome-wide contact maps. In the same way, **'G2-G2', 'G2-UA', 'UA-G2'** pairs will be used to generate the **G2** genome-wide contact maps.
All maps are then normalized using the ICE method.
