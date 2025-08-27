import argparse
from email import parser
from os import system, path, remove, makedirs, environ
import random
import string
from Bio import SeqIO
from length_bias_functions_select import *
from joblib import Parallel, delayed
import glob
from itertools import combinations, islice
import pandas as pd
import math
import hashlib
from pathlib import Path

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
                        help='number of threads to run in parallel for single-machine use or within each process_pair')
    parser.add_argument('--keep_alignments',
                        default=None,
                        action='store_true',
                        help='Whether to discard alignment files after length bias is calculated.')
    parser.add_argument('--window_size',
                        default=500,
                        type=int,
                        help='window size for diversity calculation')
    parser.add_argument('--chunks', 
                        type=int, 
                        default=None,
                        help='Total number of chunks to split pairs into (for multi-node).')
    parser.add_argument('--chunk_id', 
                        type=int, 
                        default=None,
                        help='0-based chunk id to run.')


    args = parser.parse_args()
    check_inputs(args)

    header = ['Strain 1',
                 'Strain 2',
                 'Initial divergence raw',
                 'Initial divergence',
                 'Alignment size raw',
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

    if args.chunks is not None and args.chunk_id is not None:
        genomes = sorted(glob.glob(path.join(args.genome_dir, f"*{args.genome_ext}")))
        if len(genomes) < 2:
            raise ValueError("Need ≥2 genomes")

        n = len(genomes)
        total_pairs = n * (n - 1) // 2
        chunk_id0 = args.chunk_id - 1
        start, end = contiguous_slice(total_pairs, args.chunks, chunk_id0)
        pair_iter = islice(combinations(genomes, 2), start, end)

        out_files = Parallel(n_jobs=args.num_threads,prefer="processes",batch_size=1)(
            delayed(process_pair)(g1, g2, args) for (g1, g2) in pair_iter
        )
        for f in out_files:
            print(f)
        rows = [open(f).read().strip().split() for f in out_files]
        df = pd.DataFrame(rows, columns=header)
        outfile_name = '{final_output_dir}/{base_name}.length_bias_{window_size}.txt'.format(final_output_dir=args.final_output_dir, base_name=args.base_name, window_size=args.window_size)
        df.to_csv(outfile_name, sep='\t', index=False)

        return  # <-- IMPORTANT so we don't fall through to single-node mode


    length_bias_files = run_on_single_machine(args.num_threads,
                                                  args.genome_dir,
                                                  args.genome_ext,
                                                  args.alignment_dir,
                                                  args.window_size,
                                                  args.keep_alignments,
                                                  args.final_output_dir)
    

    rows = [open(f).read().strip().split() for f in length_bias_files]
    df = pd.DataFrame(rows, columns=header)
    outfile_name = '{final_output_dir}/{base_name}.length_bias_{window_size}.txt'.format(final_output_dir=args.final_output_dir, base_name=args.base_name, window_size=args.window_size)
    df.to_csv(outfile_name, sep='\t', index=False)

def set_single_thread_env():
    # keep all native libs to 1 thread
    for var in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
                "BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        environ[var] = "1"

def check_inputs(args):

    # Check that contig files exist in the directory
    contig_list = glob.glob('{contigdir}/*{extension}'.format(contigdir=args.genome_dir,extension=args.genome_ext))

    if len(contig_list) == 0:
        raise FileNotFoundError('Files with contig extension not found in directory.')

    # Check for alignment directory. Makes it if it isn't there.
    if not path.exists(args.alignment_dir):
        print('Alignment output directory does not exist. Creating new directory.')
        makedirs(args.alignment_dir, exist_ok=True)

    # Check for final ouput_directory. Makes it if it isn't there.
    if not path.exists(args.final_output_dir):
        print('Final output directory does not exist. Creating new directory.')
        makedirs(args.final_output_dir, exist_ok=True)

def run_on_single_machine(threads,
                          genome_directory,
                          genome_extension,
                          alignment_dir,
                          window_size,
                          keep_alignments,
                          final_output_dir):
    
    renamed_genomes = [rename_for_mugsy(g) for g in glob.glob(genome_directory + '*' + genome_extension)]
    pairs_and_seeds = [(g1, g2, random.randint(1, int(1e9))) for g1, g2 in combinations(renamed_genomes, 2)]
    length_bias_files = Parallel(n_jobs=threads)(delayed(align_and_calculate_length_bias)(g1, g2, alignment_dir, window_size, seed, keep_alignments, final_output_dir) for g1, g2, seed in pairs_and_seeds)
    return length_bias_files

def contiguous_slice(n_items, chunks, chunk_id):
    """Return [start, end) for a contiguous slice of a range(n_items)."""
    if chunks is None or chunk_id is None:
        return 0, n_items
    if not (0 <= chunk_id < chunks):
        raise ValueError(f'chunk_id {chunk_id} out of range [0, {chunks-1}]')
    per = math.ceil(n_items / chunks)
    start = per * chunk_id
    end = min(n_items, start + per)
    return start, end

def deterministic_seed(g1, g2):
    """Stable per-pair seed from filenames (avoid cross-task collisions)."""
    s = (Path(g1).name + '|' + Path(g2).name).encode()
    return int.from_bytes(hashlib.sha1(s).digest()[:4], 'little')

def process_pair(g1, g2, args):
    set_single_thread_env() 
    g1r = rename_for_mugsy(g1)
    g2r = rename_for_mugsy(g2)
    seed = deterministic_seed(g1, g2)
    return align_and_calculate_length_bias(
        g1r, g2r,
        args.alignment_dir,
        args.window_size,
        seed,                    # replaces random.randint(...) for reproducibility
        args.keep_alignments,
        args.final_output_dir
    )

if __name__ == '__main__':
    main()
