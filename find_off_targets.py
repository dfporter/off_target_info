import pandas, re, collections, os, sys, argparse
from typing import List, Mapping, Union


def make_mRNA_bed_file(genome_gtf_fname: str) -> None:

    gencode = pandas.read_csv(
        genome_gtf_fname, comment="#", sep="\t",
        names=['seqname', 'source', 'feature', 'start' , 'end', 'score', 
        'strand', 'frame', 'attribute'])

    def gene_info(line):
        name, _type, status, level, tsl = ('', '', '', '', '')
        if (m := re.search('gene_name "([^"]+)"', line)) is not None:
            name = m.group(1)
        if (m := re.search('gene_type "([^"]+)"', line)) is not None:
            _type = m.group(1)
        if (m := re.search('gene_status "([^"]+)"', line)) is not None:
            status = m.group(1)
        if (m := re.search('level "([^"]+)"', line)) is not None:
            level = m.group(1)
        if (m := re.search('transcript_support_level "([^"]+)"', line)) is not None:
            tsl = m.group(1)
        return (name, _type, status, tsl)

    print(gencode.head())

    # Extract genes
    gencode_genes = gencode[(gencode.feature == "CDS")][['seqname', 'start', 'end', 'attribute']
            ].copy().reset_index().drop('index', axis=1) 

    # Set annotation columns.
    (gencode_genes["gene_name"], gencode_genes["gene_type"],
            gencode_genes["gene_status"], gencode_genes["gene_level"]) = zip(
            *gencode_genes.attribute.apply(lambda x: gene_info(x)))

    # Extract mRNA genes.
    gencode_mRNA = gencode_genes.loc[
        [x=='protein_coding' for x in gencode_genes.gene_type],:]

    # Write output in bed format.
    gencode_mRNA.loc[:,['seqname', 'start', 'end', 'gene_name']].to_csv(
        mRNA_bed_fname, header=False, index=False, sep='\t')


def write_input_to_offinder(genome_fa_fname: str, gRNA_list: List) -> None:
    
    def write_offinder_line(gRNAs: str) -> str:
        if type(gRNAs) != type(''):
            return ''
        return ''.join([f'{seq} 5\n' for n, seq in enumerate(gRNAs.split(' '))])

    offin_lines = [f'{genome_fa_fname}\nNNNNNNNNNNNNNNNNNNNNNGG\n']
    offin_lines += [write_offinder_line(gRNAs) for gRNAs in gRNA_list]
    offin_lines = ''.join([x for x in offin_lines if x!=''])

    # Write the file used by the program cas-offinder to find genomic binding sites.
    with open(cas_offinder_input_fname, 'w') as f:
        f.write(offin_lines)


def intersect(
    cas_offinder_output_fname: str,
    cas_offinder_output_as_bed_fname: str,
    mRNA_bed_fname: str,
    overlaps_bed_fname: str) -> None:
    """Write overlaps_bed_fname using mRNAs bed file and cas-offinder output.
    cas-offinder output includes the genomic location of the left border of
    the 23 nt version of the gRNA (that is, N*20 + NGG). This position is the PAM
    site if the gRNA is mapping to the reverse strand, and ~20 nt away if the gRNA
    is mapping to the positive strand."""

    gRNA_mapping_locations = pandas.read_csv(cas_offinder_output_fname, sep='\t', header=None)

    gRNA_mapping_locations.columns = ['gRNA', 'chr',  'pos', 'match', 'strand', '# mismatches']
    gRNA_mapping_locations['chr'] = [x.split(' ')[0] for x in gRNA_mapping_locations['chr']]

    # Subset mismatches:
    gRNA_mapping_locations = gRNA_mapping_locations.loc[gRNA_mapping_locations['# mismatches']<4]

    print(f'# gRNAs: {len(set(gRNA_mapping_locations.gRNA))}')
    print(gRNA_mapping_locations.shape)


    to_cut_pos = lambda pos, strand: (pos+17 if strand == '+' else pos-6)

    gRNA_mapping_locations['start'] = [
        to_cut_pos(pos, strand) for pos, strand in zip(gRNA_mapping_locations['pos'], gRNA_mapping_locations['strand'])]
    gRNA_mapping_locations['end'] = [
        to_cut_pos(pos, strand)+1 for pos, strand in zip(gRNA_mapping_locations['pos'], gRNA_mapping_locations['strand'])]

    gRNA_mapping_locations.loc[:,['chr', 'start', 'end', 'gRNA']].to_csv(
        cas_offinder_output_as_bed_fname, header=False, index=False, sep='\t')

    cmd = f'bedtools intersect -a {cas_offinder_output_as_bed_fname} -b {mRNA_bed_fname} -wb -wa'
    cmd += f' > {overlaps_bed_fname}'
    print(cmd)
    os.system(cmd)

def is_likely_target(mm_locs: List[List[int]]) -> bool:

    if len(mm_locs) == 0:
        return True
    any_good = False

    for locs in mm_locs:
        if any([x>17 for x in locs]):
            continue
        if all([x<8 for x in locs]):
            return True

    return False

def identify_target_genes(
    input_excel_df: pandas.DataFrame,
    overlaps_bed_fname: str, intended_vs_actual_targets_fname: str
    ) -> List[pandas.DataFrame]:
    """Process the output of cas-offinder (genomic locations) by 
    converting it to bed format and intersecting with RNA locations.
    """

    # Make a one-guide-per-line version of the input dataframe,
    # just holding the intended targets of each gRNA.
    info = df_of_gRNA_seq_to_targets(input_excel_df, gRNA_seq_col_name, target_gene_col_name)
    info.columns = ['Intended targets']
    print("Simplification of input file of gRNAs:\n", info.head(2))

    # Load gRNA and all actual genomic hits (output 
    # from the bedtools intersect command).
    # Positions:
    # cas-offinder_output_{sample}.bed -> hold cut position as start/end coordinates
    # cas-offinder_output_{sample}.txt -> holds left border of gRNA 23 nt mapping position.
    # overlaps_{sample}.bed -> gRNA positions are the same as the bed file (cut position).

    overlaps = pandas.read_csv(overlaps_bed_fname, sep='\t', header=None)
    overlaps.columns = ['gRNA chr', 'gRNA start', 'gRNA end', 'gRNA', 'gene chr', 
        'gene start', 'gene end', 'gene']
    overlaps.index = overlaps['gRNA'].tolist()
    overlaps.sort_values(by=['gRNA'])
    overlaps['t'] = [a for a in zip(overlaps['gRNA'], overlaps['gRNA chr'], overlaps['gRNA start'])]

    # Add mismatch locations.
    overlaps = mismatch_locations(overlaps, cas_offinder_output_fname)

    print(f"lines from cas-offinder output: {overlaps.shape[0]}")
    print(f"gRNAs in cas-offinder output: {len(set(overlaps['gRNA']))}")

    # Simplify the overlaps dataframe to just hold the gRNAs and targets.
    hits = df_of_gRNA_seq_to_targets(overlaps, 'gRNA', 'gene')
    hits['gRNA'] = list(hits.index)
    hits = pandas.DataFrame(hits.groupby('gRNA')['Targets'].apply(lambda x: set(*x)))
    print("Simplification of cas-offinder output:\n", hits.head())
    # Add the actual targets to the dataframe with intended targets.
    info['Targets'] = [hits.loc[x, 'Targets'] if x in overlaps.index else set() for x in info.index]
 
    # Transfer mismatch locations to the info dataframe.
    info.loc[:, 'mm locs'] = ''
    info.loc[:, 'Good targets'] = [set() for _ in range(len(info.index))]
    for gRNA, targets in zip(info.index, info.Targets):
        for target in list(targets):
            mm_locs = [mm for gene, mm in zip(overlaps.loc[gRNA, 'gene'], overlaps.loc[gRNA, 'mm locs']) if gene==target]
            mm_locs = list(set(tuple(_) for _ in mm_locs))  # No duplicates.
            info.loc[gRNA, 'Good targets'] |= set([target]) if is_likely_target(mm_locs) else set()
            mm_locs = re.sub('\[\[', '[', str(mm_locs))
            mm_locs = re.sub('\]\]', ']', mm_locs)
            info.loc[gRNA, 'mm locs'] += f"{target}: {mm_locs}; "

    info.to_excel(intended_vs_actual_targets_fname)

    print("With actual targets:\n", info.tail())
    print(f"Wrote output to {intended_vs_actual_targets_fname}")

    return overlaps, info


def mismatch_locations(
    overlaps: pandas.DataFrame,
    cas_offinder_output_fname: str,) -> pandas.DataFrame:

    # cas-offinder_output_{sample}.txt -> holds left border of gRNA 23 nt mapping position.
    # overlaps_{sample}.bed -> gRNA positions are the same as the bed file (cut position).
    to_cut_pos = lambda pos, strand: (pos+17 if strand == '+' else pos-6)
    
    # Process cas_offinder_output_fname into mapping locations:
    def interpret_mismatches(_str):
        def is_match(x):
            return 0 if x.isupper() else 1
        return [n for n, base in enumerate(_str) if is_match(base)]

    gRNA_mapping_locations = pandas.read_csv(cas_offinder_output_fname, sep='\t', header=None)
    gRNA_mapping_locations.columns = ['gRNA', 'chr',  'pos', 'match', 'strand', '# mismatches']

    # Need to change the position to be the cut position in order to match the overlaps
    # bed file.
    gRNA_mapping_locations['pos'] = [
        to_cut_pos(pos, strand) for pos, strand in zip(gRNA_mapping_locations['pos'], gRNA_mapping_locations['strand'])]

    gRNA_mapping_locations = gRNA_mapping_locations[gRNA_mapping_locations['# mismatches']<4]
    gRNA_mapping_locations['chr'] = [x.split(' ')[0] for x in gRNA_mapping_locations['chr']]
    gRNA_mapping_locations['mm locs'] = gRNA_mapping_locations['match'].apply(interpret_mismatches)

    index = [(str(x), str(y), z) for x,y,z in zip(gRNA_mapping_locations['gRNA'], gRNA_mapping_locations['chr'], gRNA_mapping_locations['pos'])]
    gRNA_mapping_locations.set_index(['gRNA', 'chr', 'pos'], inplace=True)

    print(gRNA_mapping_locations.head(2))
    # Finished processing cas_offinder_output_fname.

    def get_locs(gRNA_and_loc):
        if gRNA_and_loc not in gRNA_mapping_locations.index:
            print(f"Couldn't find {gRNA_and_loc}...")
            return []
        try:
            mm_locs = gRNA_mapping_locations.loc[gRNA_and_loc, 'mm locs']
        except:
            print("No 'mm locs' column in gRNA_mapping_locations file.")
            return []

        if type(mm_locs) == type([]):
            return mm_locs

        return gRNA_mapping_locations.loc[gRNA_and_loc, 'mm locs'].values[0]

    overlaps['mm locs'] = [get_locs(x) for x in overlaps['t']]

    return overlaps  # It was edited in place, but for clarity.


def df_of_gRNA_seq_to_targets(
    _df: pandas.DataFrame, gRNA_seq_col, target_col) -> pandas.DataFrame:

    q = [
        [(gRNAs.split(' ')[n], gene) for n in range(len(gRNAs.split(' ')))] for 
        gene, gRNAs in zip(_df[target_col], _df[gRNA_seq_col])]
    
    flat = collections.defaultdict(set)
    for a in q:
        [flat[gRNA_seq].add(target_name) for (gRNA_seq, target_name) in a]
        
    df = pandas.DataFrame(index=flat.keys(), columns=['Targets'])
    df['Targets'] = flat.values()

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Use cas-offinder to examine which genes are targeted by gRNAs, and the locations of mismatches.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_excel', type=str, metavar='input_excel',
        help='Excel file with gRNAs to evaluate.')
    parser.add_argument('--data_folder', type=str,
        default='data/processed/', help='Folder to hold outputs.')
    parser.add_argument('--sample', type=str,
        default='all_sources', help='Sample name for suffixing output file names.')
    parser.add_argument('--genome_gtf', type=str,
        default='/opt/genomes/gencode.v29/gencode.v29.primary_assembly.annotation.gtf', 
        help='Genomic gtf.')
    parser.add_argument('--genome_fasta', type=str,
        default="/opt/genomes/gencode.v29/GRCh38.primary_assembly.genome.fa", 
        help='Genomic fasta filename.')
    parser.add_argument('--gRNA_seq_col_name', type=str,
        default='all gRNAs',
        help='Name of the column in the input_excel file holding guide RNAs sequences (guide RNAs separated by spaces.)')
    parser.add_argument('--target_gene_col_name', type=str,
        default='Gene name',
        help='Name of the column in the input_excel file holding intended target gene names.')
    parser.add_argument('--do_not_run_cas_offinder', action='store_true',
        default=False,
        help="Don't run cas-offinder; assume its output is present. Cas-offinder takes a long time to run.")
    parser.add_argument('--overwrite_genome_bed_file', action='store_true',
        default=False,
        help="Overwrite a genome bed file if present.")
    parser.add_argument('--cas_offinder_path', 
        default='cas-offinder',
        help="Path to the cas-offinder binary executable.")
    
    args = parser.parse_args()


    # File paths.
    data_folder = args.data_folder
    sample = args.sample
    genome_gtf_fname = args.genome_gtf
    genome_fa_fname = args.genome_fasta
    input_excel = args.input_excel
    gRNA_seq_col_name = args.gRNA_seq_col_name
    target_gene_col_name = args.target_gene_col_name
    cas_offinder_path = args.cas_offinder_path

    # Make the datafolder if necessary.
    os.makedirs(data_folder, exist_ok=True)

    # Settings.
    run_cas_offinder = not(args.do_not_run_cas_offinder)
    overwrite_genome_bed_file = args.overwrite_genome_bed_file

    # Generated file paths for intermediates.
    cas_offinder_input_fname = f'{data_folder}/cas-offinder_input_{sample}.txt'
    cas_offinder_output_fname = f'{data_folder}/cas-offinder_output_{sample}.txt'
    cas_offinder_output_as_bed_fname = f'{data_folder}/cas-offinder_output_{sample}.bed'
    mRNA_bed_fname = f'{data_folder}/mRNAs.bed'
    overlaps_bed_fname = f'{data_folder}/overlaps_{sample}.bed'

    # Output file.
    intended_vs_actual_targets_fname = f'{data_folder}/intended_and_actual_targets_{sample}.xlsx'

    input_excel_df = pandas.read_excel(input_excel)
    input_excel_df = input_excel_df.loc[
        [type(n)==type('') for n in input_excel_df[gRNA_seq_col_name]]]
    write_input_to_offinder(genome_fa_fname, input_excel_df[gRNA_seq_col_name])

    # Write the genome bed file for all mRNAs.
    if (not os.path.exists(mRNA_bed_fname)) or overwrite_genome_bed_file:
        make_mRNA_bed_file(genome_gtf_fname)

    # Running cas-offinder.
    # C denotes the use of CPUs vs GPUs. This program takes >=20 minutes.
    cmd = f"{cas_offinder_path} {cas_offinder_input_fname} C {cas_offinder_output_fname}"
    print(cmd)
    run_cas_offinder and os.system(cmd)

    # Assign mapping locations to genes.
    intersect(
        cas_offinder_output_fname, cas_offinder_output_as_bed_fname,
        mRNA_bed_fname, overlaps_bed_fname)

    # Write the final output.
    identify_target_genes(
        input_excel_df, overlaps_bed_fname, intended_vs_actual_targets_fname)

