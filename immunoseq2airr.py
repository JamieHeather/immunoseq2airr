# -*- coding: utf-8 -*-

"""
immunoseq2airr.py

Convert v2-exported data from the immunoSEQ/immuneACCESS platforms (Adaptive Biotech)
into the AIRR Community rearrangement format proposed in Vander Heiden et al, 2018, Frontiers in Immunology

"""

import gzip
import argparse
import os

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '1.0.0'
__author__ = 'Jamie Heather'


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
      description="immunoseq2airr v" + str(__version__) +
                  ": Convert v2 Adaptive TCR data exported data into the AIRR Community rearrangement format")

    # Input and output options
    parser.add_argument('-i', '--input', required=True, type=str,
                        help='Path to input file. Required.')

    parser.add_argument('-n', '--name', required=False, type=str,
                        help='String to use as prefix in sequence_id field. Optional; default is input file name.')

    parser.add_argument('-o', '--output', required=False, type=str,
                        help='Path to output file. Optional; if not specified will use modified input filename.')

    parser.add_argument('-lz', '--leading_zero', required=False, type=int, default=8,
                        help='How many digits to fill leading zeroes in sequence_id field do. Optional. Default = 8.')

    parser.add_argument('-z', '--compress', action='store_true', required=False,
                        help='Option to compress output file.')

    parser.add_argument('-or', '--ignore_orphons', action='store_true', required=False,
                        help='Ignore orphons (genes on other chromosomes, e.g. TRBV20/OR9-2) in ambiguous gene calls')

    parser.add_argument('-d', '--no_d_ambiguity', action='store_true', required=False,
                        help='Prevent output of ambiguous gene calls, e.g. \'TRBD1,TRBD2\' would just be ignored')

    parser.add_argument('-nd', '--no_d', action='store_true', required=False,
                        help='Ignore all D information. Enabled by default if processing alpha chain.')

    parser.add_argument('-c', '--chain', required=False, type=str, default='b',
                        help='TCR chain. Alpha/A/TRA/a or Beta/B/TRB/b. Default = b.')

    parser.add_argument('-a', '--allow_ambiguity', action='store_true', required=False,
                        help='Allow ambiguity where V/J allele information not given. Just output whatever\'s there.' +
                        '\n\tNB: check output if using this option, as it may mask data in an unexpected format.')

    parser.add_argument('-af', '--abundance_filter', required=False, type=float,
                        help='Optional abundance filter. Only retain rearrangements with a value equal or greater than')

    parser.add_argument('-pf', '--productivity_filter', action='store_true', required=False,
                        help='Optional productivity filter. Only retain potentially productive rearrangements')

    parser.add_argument('-sa', '--strip_alleles', action='store_true', required=False,
                        help='Optionally strip allele information from gene calls (e.g. TRBV5 not TRBV5*01)')

    parser.add_argument('-mf', '--motif_filter', action='store_true', required=False,
                        help='Optional motif filter. Discard TCRs with CDR3s not bound by C and F (or W for TRA)')

    # TODO add option to specify output directory?

    return parser.parse_args()


def opener(in_file):
    """
    Gets the appropriate opener for an input file, and checks it exists
    :param in_file: Path to a v2-formatted immunoSEQ/immuneACCESS file
    :return: The file opened with the appropriate opening command
    """
    # First check file exists
    if os.path.isfile(in_file):

        # If it does, open it and return the opened file
        if in_file.endswith('.gz'):
            return gzip.open(in_file)
        elif in_file.endswith('.tsv'):
            return open(in_file)
        else:
            raise IOError("Error, unknown file extension: " + in_file)

    else:
        raise IOError("Error, cannot locate file: " + in_file)


def get_name(in_file):
    """
    :param in_file: Path to file to convert
    :return: The inferred sample name, defined by file name shorn of any file extensions
    """

    return in_file.split('/')[-1].split('.')[0]


def tidy_gene(tcr_gene):
    """
    :param tcr_gene: The named TCR gene as read in the data file
    :return: A properly IMGT-recognised TCR gene name
    """

    # Remove leading zeroes and inappropriate 'C' in 'TRXY'
    tcr_gene = tcr_gene.replace('TCR', 'TR').replace('V0', 'V').replace('J0', 'J').replace('D0', 'D').replace('-0', '-')

    # Correct orphon names
    tcr_gene = tcr_gene.replace('-or09_02', '/OR9-2')

    # Remove inappropriate sub-family naming in TRAJ and TRBD genes
    if 'TRAJ' in tcr_gene or 'TRBD' in tcr_gene:
        tcr_gene = tcr_gene.replace('-1', '')

    return tcr_gene


def check_gene(gene_call):
    """
    Check to make sure a novel Adaptive gene name gets used
    :param gene_call: the gene name as reported (NB, could either have or lack allele level information)
    :return: the gene name as proscribed by IMGT (with allele level info if present)
    """

    gene_bits = gene_call.split('*')

    if len(gene_bits) == 1:
        if gene_call in adaptive_v_convert:
            return adaptive_v_convert[gene_call]
        else:
            return gene_call

    elif len(gene_bits) == 2:
        if gene_bits[0] in adaptive_v_convert:
            return adaptive_v_convert[gene_bits[0]] + '*' + gene_bits[1]
        else:
            return gene_call

    else:
        raise IOError("Inappropriate gene format detected for gene name - " + gene_call)


def convert_tcr(split_line, tcr_id):
    """
    Performs the actual conversions
    :param split_line: the line read in from the input v2 tsv, appropriately opened/rstripped/split on tabs
    :param tcr_id: a unique identifier for a given rearrangement, in this case a zero-padded number prefixed by a name
    :return:
    """

    # Compile the necessary fields for output
    out_vals = {'sequence_id': tcr_id,
                'sequence': split_line[params['sequence_index']], 'rev_comp': 'F',
                'duplicate_count': split_line[params['abundance_index']]}

    # If the option has been set, only retain those sequences with a value equal to or greater than that threshold
    if input_args['abundance_filter']:
        if float(out_vals['duplicate_count']) < input_args['abundance_filter']:
            return

    # Infer productivity (using presence of CDR3 and Adaptive sequenceStatus value), take junction if there
    if split_line[params['cdr3_index']] and split_line[params['productivity']] == 'In':
        out_vals['junction_aa'] = split_line[params['cdr3_index']]
        out_vals['productive'] = 'T'

        # If the option to discard rearrangements lacking proper CDR3 motifs has been set, skip this entry if not C/F
        if input_args['motif_filter']:
            if out_vals['junction_aa'][0] != 'C':
                return
            elif chain == 'TRB' and out_vals['junction_aa'][-1] != 'F':
                return
            # Human TRAJ are a bit more flexible as to their junction-defining residue
            elif chain == 'TRA' and out_vals['junction_aa'][-1] not in ['F', 'W', 'C']:
                return

    else:
        out_vals['junction_aa'] = ''
        out_vals['productive'] = 'F'

        # If the option to ignore non-productive rearrangements has been set, skip this row
        if input_args['productivity_filter']:
            return

    # TODO add a CDR3 sanity length/motif filter?

    # If users wanted to they could infer the junction nt sequence, but I haven't, as it's redundant/not very useful
    out_vals['junction'] = ''

    # Extract the VDJ genes, fixing their nomenclature and combining together multiple possible calls
    for gene in ['v', 'd', 'j']:

        # First take Adaptive's best call
        call = split_line[params[gene + 'MaxResolved']]

        # Check whether the code wants to be looking for D genes
        if input_args['no_d'] and gene == 'd':
            sorted_call = ''

        # Check whether a gene has been called - if not (and not D) check in the ambiguous gene name ties field
        # NB ambiguous Ds are ignored by default, as there are only two options for TRBD and they're almost identical
        elif not call or call == 'unresolved':
            if gene == 'd' and input_args['no_d_ambiguity']:
                sorted_call = ''
            else:
                sorted_call = resolve_ambiguous_name(split_line, gene)

        # If it has full allele accuracy (indicated by an asterisk), tidy it up and take that as the result
        elif call[-3] == '*':
            sorted_call = check_gene(tidy_gene(call))

        # Depending on the (hidden) version of the input data, remaining ambiguity might be resolved in 2 places:
        # either in the GeneNameTies or AlleleNameTies fields - need to infer which is correct and deal appropriately
        else:

            if bits[params[gene + 'GeneNameTies']]:
                sorted_call = resolve_ambiguous_name(split_line, gene)
            elif bits[params[gene + 'GeneAlleleTies']]:
                sorted_call = resolve_ambiguous_allele(call, split_line, gene)

            # However some files are not even covered by that broad formatting, so you just need to allow whatever
            elif input_args['allow_ambiguity']:
                sorted_call = check_gene(tidy_gene(call))

            else:
                raise IOError("Unknown format on line " + str(line_count) + "! Cannot continue. "
                     "\n\tAmbiguity for " + gene.upper() + " gene calls lacking allele info that is"
                     "\n\t not resolved in either 'Gene' or 'Allele Ties' fields."
                     "\n\tTry re-running the script using the '-a' flag (to allow ambiguity),"
                     "\n\t and check that the format of the output document is correct.")

        # If option is selected, remove allele level information
        if input_args['strip_alleles']:
            sorted_call = strip_alleles(poss_alleles, sorted_call)

        out_vals[gene + '_call'] = sorted_call

    # Finally pad the missing values for the required columns
    for value in [x for x in out_headers if x not in out_vals]:
        out_vals[value] = ''

    return out_vals


def resolve_ambiguous_name(line_split, current_gene):
    """
    Resolve ambiguous gene calls where there is multiple gene level information available
    NB requires global 'params' dict, and global 'input_args' dict (for 'ignore_orphons'] declaration)
    :param line_split: 'bits', i.e. the rstripped line of the tsv split on tab characters
    :param current_gene: v/d/j
    :return: the sorted call (all possible gene calls joined by commas)
    """

    unresolved = line_split[params[current_gene + 'GeneNameTies']].split(',')

    # If users have selected, ignore orphons when compiling ambiguous gene possibilities
    if input_args['ignore_orphons']:
        unresolved = [x for x in unresolved if 'or' not in x]

    return','.join([check_gene(tidy_gene(x)) for x in unresolved])


def resolve_ambiguous_allele(current_call, line_split, current_gene):
    """
    Resolve ambiguous gene calls where there is multiple allele level information available, for one gene
    NB requires global 'params' dict
    :param current_call: What the gene is 'currently' determined to be
    :param line_split: 'bits', i.e. the rstripped line of the tsv split on tab characters
    :param current_gene: v/d/j
    :return: the sorted call (all possible gene calls joined by commas)
    """

    return ','.join([check_gene(tidy_gene(current_call)) + '*' + x
                     for x in line_split[params[current_gene + 'GeneAlleleTies']].split(',')])


def infer_chain(chain_str):
    """
    :param chain_str: An input string, used to refer to a specific TCR chain locus
    :return: TRA or TRB
    """

    if chain_str.upper() in ['TRA', 'A', 'ALPHA']:
        return 'TRA'
    elif chain_str.upper() in ['TRB', 'B', 'BETA']:
        return 'TRB'
    else:
        raise IOError("Chain string not recognised: please use TRA/B, A/B, or alpha/beta")


def strip_alleles(possible_alleles, gene_call):
    """
    :param possible_alleles: list of possible allele strings (e.g. '*01, *02'...). Default is to *09, more than in IMGT
    :param gene_call: The current working call for the gene(s) in a given rearrangement
    :return: The same gene(s) stripped of any allelic information
    """

    # First remove allele strings
    for pa in possible_alleles:
        gene_call = gene_call.replace(pa, '')

    # Then ensure we're still dealing with unique gene calls
    genes = list(set(gene_call.split(',')))
    genes.sort()

    return ','.join(genes)


# Need to account for the fact that Adaptive have taken it on themselves to rename the TCR genes
# TODO add mappings for other loci to make applicable to all chains *******
# TODO transfer to supplementary file to make it easier to change?
adaptive_v_convert = {'TRBV1-1': 'TRBV1',
                      'TRBV2-1': 'TRBV2',
                      'TRBV9-1': 'TRBV9',
                      'TRBV13-1': 'TRBV13',
                      'TRBV14-1': 'TRBV14',
                      'TRBV15-1': 'TRBV15',
                      'TRBV16-1': 'TRBV16',
                      'TRBV17-1': 'TRBV17',
                      'TRBV18-1': 'TRBV18',
                      'TRBV19-1': 'TRBV19',
                      'TRBV26-1': 'TRBV26',
                      'TRBV27-1': 'TRBV27',
                      'TRBV28-1': 'TRBV28',
                      'TRBV30-1': 'TRBV30',
                      'TRAV14-1': 'TRAV14/DV4',
                      'TRAV23-1': 'TRAV23/DV6',
                      'TRAV29-1': 'TRAV29/DV5',
                      'TRAV36-1': 'TRAV36/DV7',
                      'TRAV10-1': 'TRAV10',
                      'TRAV11-1': 'TRAV11',
                      'TRAV15-1': 'TRAV15',
                      'TRAV16-1': 'TRAV16',
                      'TRAV17-1': 'TRAV17',
                      'TRAV18-1': 'TRAV18',
                      'TRAV19-1': 'TRAV19',
                      'TRAV2-1': 'TRAV2',
                      'TRAV20-1': 'TRAV20',
                      'TRAV21-1': 'TRAV21',
                      'TRAV22-1': 'TRAV22',
                      'TRAV24-1': 'TRAV24',
                      'TRAV25-1': 'TRAV25',
                      'TRAV27-1': 'TRAV27',
                      'TRAV28-1': 'TRAV28',
                      'TRAV3-1': 'TRAV3',
                      'TRAV30-1': 'TRAV30',
                      'TRAV31-1': 'TRAV31',
                      'TRAV32-1': 'TRAV32',
                      'TRAV33-1': 'TRAV33',
                      'TRAV34-1': 'TRAV34',
                      'TRAV35-1': 'TRAV35',
                      'TRAV37-1': 'TRAV37',
                      'TRAV39-1': 'TRAV39',
                      'TRAV4-1': 'TRAV4',
                      'TRAV40-1': 'TRAV40',
                      'TRAV41-1': 'TRAV41',
                      'TRAV46-1': 'TRAV46',
                      'TRAV5-1': 'TRAV5',
                      'TRAV6-1': 'TRAV6',
                      'TRAV7-1': 'TRAV7',
                      'TRDV1-1': 'TRDV1',
                      'TRDV2-1': 'TRDV2',
                      'TRDV3-1': 'TRDV3'
                      }

# TODO potentially infer parameter columns from headers (or even early rows)?
# Would help with rogue file versions, but would probably require knowing all the file formats in the first place
params = {'sequence_index': 0, 'cdr3_index': 1, 'abundance_index': 2, 'productivity': 38,
          'vMaxResolved': 5, 'dMaxResolved': 12, 'jMaxResolved': 19,
          # 'vFamilyName': 6, 'dFamilyName': 13, 'jFamilyName': 20,
          'vGeneNameTies': 10, 'dGeneNameTies': 17, 'jGeneNameTies': 24,
          'vGeneAlleleTies': 11, 'dGeneAlleleTies': 18, 'jGeneAlleleTies': 25}

# AIRR-seq community format out headers
out_headers = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'junction_aa', 'duplicate_count',
               'rev_comp', 'productive', 'sequence_alignment', 'germline_alignment',
               'junction', 'v_cigar', 'd_cigar', 'j_cigar']

if __name__ == '__main__':

    # Get input command line arguments
    input_args = vars(args())
    chain = infer_chain(input_args['chain'])
    if chain == 'TRA':
        input_args['no_d'] = True

    # Get names to use for output rearrangements/file
    name = get_name(input_args['input'])
    if input_args['name']:
        name = input_args['name']

    if input_args['output']:
        output_name = input_args['output'] + '.tsv'
    else:
        output_name = name + '.tsv'
        if output_name == input_args['input'].split('/')[-1]:
            output_name = name + '-airr.tsv'

    if input_args['strip_alleles']:
        # Establish which allele strings to remove from gene call strings if requested
        # There are currently only up to *07 stored in IMGT, so this is a generous allowance of possible alleles
        poss_alleles = ['*' + str(x).zfill(2) for x in range(1, 9)]

    # Check for need to compress  # TODO ideally switch to a splittable compression format
    if input_args['compress']:
        output_name += '.gz'
        out_opener = gzip.open
    else:
        out_opener = open

    # Then run through input file and process
    with out_opener(output_name, 'w') as out_file:
        line_count = 0
        out_file.write('\t'.join(out_headers) + '\n')
        out_str = ''

        for line in opener(input_args['input']):
            bits = line.rstrip().split('\t')

            # Take headers from very first line
            if line_count == 0:
                headers = bits
                line_count += 1
                continue

            # Sanity check on one of the right-most columns for the first data row
            elif line_count == 1:
                if bits[params['productivity']] not in ['In', 'Out', 'Stop']:
                    raise IOError("Unexpected value in \'sequenceStatus\' column: " + bits[params['productivity']])

            # Then proceed to convert the data into the AIRR-seq community format!
            rearrangement_name = name + '|' + str(line_count).zfill(input_args['leading_zero'])
            converted = convert_tcr(bits, rearrangement_name)
            if converted:
                out_str += '\t'.join([converted[x] for x in out_headers]) + '\n'

            line_count += 1

        out_file.write(out_str)
