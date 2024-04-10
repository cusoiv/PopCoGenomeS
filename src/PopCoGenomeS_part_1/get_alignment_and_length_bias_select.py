import argparse
from os import system, path, remove, makedirs
import random
import string
from Bio import SeqIO
from length_bias_functions_select import *
from joblib import Parallel, delayed
import glob
from itertools import combinations
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=('Align contigs in a job array'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )


    parser.add_argument('--genome_dir',
                        default=None,
                        type=str,
                        help='Directory containing genome files.')
    parser.add_argument('--genome_ext',
                        default=None,
                        type=str,
                        help='Extension for genome files (e.g., .fasta')
    parser.add_argument('--alignment_dir',
                        default=None,
                        type=str,
                        help='Directory for alignments. Please provide absolute path.')
    parser.add_argument('--base_name',
                        default='output',
                        type=str,
                        help='base output file name')
    parser.add_argument('--final_output_dir',
                        default='./',
                        type=str,
                        help='Directory for final output.')
    parser.add_argument('--num_threads',
                        default=1,
                        type=int,
                        help='number of threads to run in parallel for single-machine use (i.e., not slurm)')
    parser.add_argument('--keep_alignments',
                        default=None,
                        action='store_true',
                        help='Whether to discard alignment files after length bias is calculated.')
    parser.add_argument('--window_size',
                        default=500,
                        type=int,
                        help='window size for diversity calculation')
    parser.add_argument('--hmm_cutoff',
                        default=100,
                        type=float,
                        help='cutoff to run hmm')

    args = parser.parse_args()
    check_inputs(args)


    length_bias_files = run_on_single_machine(args.num_threads,
                                                  args.genome_dir,
                                                  args.genome_ext,
                                                  args.alignment_dir,
                                                  args.window_size,
                                                  args.hmm_cutoff,
                                                  args.keep_alignments,
                                                  args.final_output_dir)
    header = ['Strain 1',
                 'Strain 2',
                 'Initial divergence iter1',
                 'Initial divergence',
                 'Alignment size iter1',
                 'Alignment size',
                 'Genome 1 size',
                 'Genome 2 size',
                 'Observed SSD',
                 'SSD 95 CI low',
                 'SSD 95 CI high',
                 'mu_div',
                 'mnb_div',
                 'sim_fr',
                 'nb_size',
                 'hmm_fr',
                 'hmm_mu',
                 'hmm_mnb']

    rows = [open(f).read().strip().split() for f in length_bias_files]
    df = pd.DataFrame(rows, columns=header)
    outfile_name = '{final_output_dir}/{base_name}.length_bias_{window_size}.txt'.format(final_output_dir=args.final_output_dir, base_name=args.base_name, window_size=args.window_size)
    df.to_csv(outfile_name, sep='\t', index=False)

def check_inputs(args):

    # Check that contig files exist in the directory
    contig_list = glob.glob('{contigdir}/*{extension}'.format(contigdir=args.genome_dir,extension=args.genome_ext))

    if len(contig_list) == 0:
        raise FileNotFoundError('Files with contig extension not found in directory.')

    # Check for alignment directory. Makes it if it isn't there.
    if not path.exists(args.alignment_dir):
        print('Alignment output directory does not exist. Creating new directory.')
        makedirs(args.alignment_dir)

    # Check for final ouput_directory. Makes it if it isn't there.
    if not path.exists(args.alignment_dir):
        print('Final output directory does not exist. Creating new directory.')
        makedirs(args.final_output_dir)

def run_on_single_machine(threads,
                          genome_directory,
                          genome_extension,
                          alignment_dir,
                          window_size,
                          hmm_cutoff,
                          keep_alignments,
                          final_output_dir):
    
    renamed_genomes = [rename_for_mugsy(g) for g in glob.glob(genome_directory + '*' + genome_extension)]
    pairs_and_seeds = [(g1, g2, random.randint(1, int(1e9))) for g1, g2 in combinations(renamed_genomes, 2)]
    length_bias_files = Parallel(n_jobs=threads)(delayed(align_and_calculate_length_bias)(g1, g2, alignment_dir, window_size, hmm_cutoff, seed, keep_alignments, final_output_dir) for g1, g2, seed in pairs_and_seeds)
    return length_bias_files

if __name__ == '__main__':
    main()
