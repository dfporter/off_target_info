# off_target_info
Use cas-offinder to examine which genes are targeted by gRNAs, and the locations of mismatches.

```txt
usage: find_off_targets.py [-h] [--data_folder DATA_FOLDER] [--sample SAMPLE] [--genome_gtf GENOME_GTF]
                           [--genome_fasta GENOME_FASTA] [--gRNA_seq_col_name GRNA_SEQ_COL_NAME]
                           [--target_gene_col_name TARGET_GENE_COL_NAME] [--do_not_run_cas_offinder]
                           [--overwrite_genome_bed_file]
                           input_excel

Use cas-offinder to examine which genes are targeted by gRNAs, and the locations of mismatches.

positional arguments:
  input_excel           Excel file with gRNAs to evaluate.

optional arguments:
  -h, --help            show this help message and exit
  --data_folder DATA_FOLDER
                        Folder to hold outputs. (default: data/processed/)
  --sample SAMPLE       Sample name for suffixing output file names. (default: all_sources)
  --genome_gtf GENOME_GTF
                        Genomic gtf. (default:
                        /opt/genomes/gencode.v29/gencode.v29.primary_assembly.annotation.gtf)
  --genome_fasta GENOME_FASTA
                        Genomic fasta filename. (default:
                        /opt/genomes/gencode.v29/GRCh38.primary_assembly.genome.fa)
  --gRNA_seq_col_name GRNA_SEQ_COL_NAME
                        Name of the column in the input_excel file holding guide RNAs sequences (guide RNAs
                        separated by spaces.) (default: all gRNAs)
  --target_gene_col_name TARGET_GENE_COL_NAME
                        Name of the column in the input_excel file holding intended target gene names. (default:
                        Gene name)
  --do_not_run_cas_offinder
                        Don't run cas-offinder; assume its output is present. Cas-offinder takes a long time to
                        run. (default: False)
  --overwrite_genome_bed_file
                        Overwrite a genome bed file if present. (default: False)

```



Example usage:
```bash
python find_off_targets.py ../Histones_pared_down.xlsx 
```
