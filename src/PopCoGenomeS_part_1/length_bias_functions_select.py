import numpy as np
from collections import Counter
from os import system, path, remove, makedirs
import random
import string
from Bio import SeqIO
from itertools import combinations, groupby
import copy
from scipy.stats import binom
from scipy.stats import poisson
from scipy.stats import nbinom
import scipy.optimize as opt
from pomegranate import *
import pandas as pd


def align_and_calculate_length_bias(genome_1_file,
                                    genome_2_file,
                                    alignment_dir,
                                    window_size, 
                                    hmm_cutoff,
                                    random_seed,
                                    keep_alignments,
                                    final_output_dir):

    alignment_file = align_genomes(genome_1_file,
                                   genome_2_file,
                                   alignment_dir,
                                   random_seed)

    length_bias_file = alignment_file + '.length_bias.txt'
    #block_dg_file_1 = alignment_file + '.block_dg_1.txt'
    calculate_length_bias(alignment_file,
                          genome_1_file,
                          genome_2_file,
                          window_size,
                          hmm_cutoff,
                          length_bias_file,
   #                       block_dg_file_1,
                          final_output_dir)
    
    if not keep_alignments:
        remove(alignment_file)
    return length_bias_file


def rename_for_mugsy(genome):
    # Assumes the strain name is everything except the extension
    strain_name = '.'.join(path.basename(genome).split('.')[0:-1])

    # We want to remove all periods and colons from sequence input so that mugsy doesn't break
    mugsy_outname = genome + '.renamed.mugsy'
    
    # Removes all bad characters
    mugsy_name = strain_name.translate(({ord(c): '_' for c in """ !@#$%^&*()[]{};:,./<>?\|`"'~-=+"""}))
    
    mugsy_s = []
    for i, s in enumerate(SeqIO.parse(genome, 'fasta')):
        s.description = ''
        s.id = '{id}_{contig_num}'.format(id=mugsy_name, contig_num=str(i))
        mugsy_s.append(s)
    SeqIO.write(mugsy_s, mugsy_outname, 'fasta')
    return mugsy_outname


def align_genomes(contig1,
                  contig2,
                  alignment_dir,
                  random_seed):
    
    random.seed(random_seed)
    # Assumes that files are named strain.contigextension.renamed.mugsy
    strain1 = '.'.join(path.basename(contig1).split('.')[0:-3])
    strain2 = '.'.join(path.basename(contig2).split('.')[0:-3])
    correct_name = '{strain1}_@_{strain2}.maf'.format(strain1 = strain1, strain2 = strain2) 
    final_name = alignment_dir+'/'+correct_name

    if not path.exists(final_name): # Only run the alignment if the file doesn't exist
        # make a temporary contig file due to parallelization issues with reading from the same filename
        out_id_1 = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))
        out_id_2 = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))


        system('cp {contig1} {alignment_dir}/{randomcontigname1}.tempcontig'.format(contig1=contig1, randomcontigname1=out_id_1, alignment_dir=alignment_dir))
        system('cp {contig2} {alignment_dir}/{randomcontigname2}.tempcontig'.format(contig2=contig2, randomcontigname2=out_id_2, alignment_dir=alignment_dir))

        # Aligning the genomes
        prefix = out_id_1 + out_id_2
        print('mugsy --directory {align_directory} --prefix {prefix} {align_directory}/{randomcontigname1}.tempcontig {align_directory}/{randomcontigname2}.tempcontig'.format(
                                                                                                                            align_directory=alignment_dir,
                                                                                                                            prefix = prefix,
                                                                                                                            randomcontigname1=out_id_1,
                                                                                                                            randomcontigname2 = out_id_2))
        system('mugsy --directory {align_directory} --prefix {prefix} {align_directory}/{randomcontigname1}.tempcontig {align_directory}/{randomcontigname2}.tempcontig'.format(
                                                                                                                            align_directory=alignment_dir,
                                                                                                                            prefix = prefix,
                                                                                                                            randomcontigname1=out_id_1,
                                                                                                                            randomcontigname2 = out_id_2))


        # Remove unneeded files
        remove('{align_directory}/{random_contig1}.tempcontig'.format(random_contig1=out_id_1, align_directory=alignment_dir))
        remove('{align_directory}/{random_contig2}.tempcontig'.format(random_contig2=out_id_2, align_directory=alignment_dir))
        if path.exists('{align_directory}/{prefix}'.format(prefix=prefix, align_directory=alignment_dir)):       #HERBOLDADDITION
            remove('{align_directory}/{prefix}'.format(prefix=prefix, align_directory=alignment_dir))
        remove('{prefix}.mugsy.log'.format(prefix=prefix))

        system('mv {random_alignment_name} {correct_name}'.format(random_alignment_name=alignment_dir+'/'+prefix +'.maf',
                                                                  correct_name=alignment_dir+'/'+correct_name))
        
    return final_name

def calculate_length_bias(input_alignment,
                          genome_1_file,
                          genome_2_file,
                          window_size,
                          hmm_cutoff,
                          output_file_1,
#                          output_file_2,
                          final_output_dir):

    g1size = sum([len(s) for s in SeqIO.parse(genome_1_file, 'fasta')])
    g2size = sum([len(s) for s in SeqIO.parse(genome_2_file, 'fasta')])

    if not path.exists(output_file_1): # only calculate the length bias if the file doesn't exist
        edge, block_dg_1 = get_transfer_measurement(input_alignment,
                             g1size=g1size,
                             g2size=g2size,
                             window_size=window_size,
                             hmm_cutoff=hmm_cutoff,
                             final_output_dir=final_output_dir)
        with open(output_file_1, 'w') as outfile:
            outfile.write(edge + '\n')
  #      with open(output_file_2, 'w') as outfile:
  #          outfile.write(block_dg_1 + '\n')


def get_transfer_measurement(alignment,
                             g1size,
                             g2size,
                             window_size,
                             hmm_cutoff,
                             final_output_dir,
                             min_block_size=0,      #HERBOLD CHANGE
#                             filtering_window=1000):
                             filtering_SNP=1):

    # Initializes local variables
    filtered_blocks = []
    strain1, strain2 = alignment.split('/')[-1].split('_@_')
    strain2 = strain2.replace('.maf', '')
    all_blocks, prefilter_total_len = get_concatenated_alignment(alignment)
    window_size=int(window_size)
    hmm_cutoff=float(hmm_cutoff)

    # Filter alignment to split into subblocks at any point where there are at least 2 gaps
    for prefilter_s1, prefilter_s2 in all_blocks:
        filtered_blocks += filter_block(prefilter_s1, prefilter_s2)
    filtered_blocks = [block for block in filtered_blocks if len(block[0]) > min_block_size]
    if len(filtered_blocks)==0:                                                                  #HERBOLD ADDITION
        edge = '\t'.join([strain1,                                                               #HERBOLD ADDITION
                        strain2,                                                                 #HERBOLD ADDITION
                        str(1),                                                                  #HERBOLD ADDITION
                        str(0),                                                                  #HERBOLD ADDITION
                        str(g1size),                                                             #HERBOLD ADDITION
                        str(g2size),                                                             #HERBOLD ADDITION
                        str(0),                                                                  #HERBOLD ADDITION
                        str(0),                                                                  #HERBOLD ADDITION
                        str(0)])                                                                 #HERBOLD ADDITION
        lb = '\t'.join([strain1,
                        strain2,
                        str(0),
                        str(0)])

    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = '1'.join(s1temp)
    Concat_S2 = '0'.join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = init_div_count * 1.0 / alignment_size

    init_div_1 = copy.copy(init_div)
    alignment_size_1 = copy.copy(alignment_size)
    concat_S1S2=[strain1,strain2,Concat_S1,Concat_S2]

    #calculate distribution of number of snps before any removal of snps
    block_div_genome_1 = block_div_calc(Concat_S1, Concat_S2, window_size)
    block_div_genome_str = ','.join(str(v) for v in block_div_genome_1)
    block_dg_1 = '\t'.join([strain1,
                          strain2,
                          block_div_genome_str])

    
    #filtering by divergence until shrinking becomes minor 
    alignment_size_old = alignment_size*2
    epsilon=1e-6
    iter=1
    
    while alignment_size_old-alignment_size>alignment_size_old*epsilon and init_div>0:
        init_div_old = copy.copy(init_div)
        alignment_size_old = copy.copy(alignment_size)
    
        final_filtered = []
    
        winlen_temp=max(int(filtering_SNP//init_div_old),100)
        winlen=min(winlen_temp,10000)
    
        slide=50
        winnumber=sum((len(block[0])//slide+1) for block in filtered_blocks)
        binom_cutoff = binom.ppf(q=1-0.05/winnumber,n=winlen,p=init_div,loc=0)
#        print (binom_cutoff)
        for s1, s2 in filtered_blocks:
            final_filtered += filter_block_by_divergence(s1, s2, init_div, winlen, winnumber, binom_cutoff,slide)
        filtered_blocks = [block for block in final_filtered if len(block[0]) > min_block_size]
        s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
        Concat_S1 = '1'.join(s1temp)
        Concat_S2 = '0'.join(s2temp)
        alignment_size = len(Concat_S1)
        init_div_count = naive_div_count(Concat_S1, Concat_S2)
        init_div = (init_div_count * 1.0) / alignment_size
#        print (init_div)
#        print (alignment_size)
        iter=iter+1
#        if iter==3:
#            break

    
    # Again scan across the concatenated genome in defined size window blocks for divergence
    block_div_genome = block_div_calc(Concat_S1, Concat_S2, window_size)
    
    if init_div > 0:
        initial = id_var_window_counts(Concat_S1, Concat_S2)
        initial_cumulative = get_cumulative_window_spectrum(initial, alignment_size)
        null_expect = single_param_null_model(np.arange(0, len(initial_cumulative)), init_div)
        observed_sum_sq_diff = np.sum(np.square(np.subtract(initial_cumulative, null_expect)))

        # Given a distribution of identical windows, bootstraps to find
        # length bias (SSD) confidence interval
        ssds = []
        for t in range(0, 200):
            initial_boot = np.random.choice(initial, len(initial), replace=True)
            initial_cumulative_boot = get_cumulative_window_spectrum(initial_boot, alignment_size)
            ssd_boot = np.sum(np.square(np.subtract(initial_cumulative_boot, null_expect)))
            ssds.append(ssd_boot)
        low_percentile = np.percentile(ssds, 0.5)
        high_percentile = np.percentile(ssds, 99.5)
    else:
        observed_sum_sq_diff = np.nan
        low_percentile = np.nan
        high_percentile = np.nan
#    initial_cumulative_str = ','.join(str(v) for v in initial_cumulative.tolist()[1:100000])
#    null_expect_str = ','.join(str(v) for v in null_expect.tolist()[1:100000])
    
    #calculate fraction of recombination using MLE, if alignment fraction >30%
    gsmall_size=min(g1size,g2size)
    alignment_fraction=alignment_size_1/gsmall_size
    if alignment_fraction>=0.3:
        mu,mnb,a,fc=P_mle(block_div_genome_1,block_div_genome,window_size)
    else:
        observed_sum_sq_diff=0
        low_percentile=0
        high_percentile=0
        fc=100
        mu=0
        mnb=0
        a=1
    mle_output=[mu,mnb,a,fc]
 #   print ([strain1,strain2])
 #   print (mle_output)
    location_dict={}
    hmm_fr=np.nan
    hmm_mu=np.nan
    hmm_mnb=np.nan
    if mu/window_size < hmm_cutoff and mu!=0:
        RC=run_hmm(block_div_genome_1,window_size,mle_output)
        RC.pop(0)
        hmm_fc=RC.count('C')/len(RC)*100
        hmm_fr=RC.count('R')/len(RC)*100
        location_dict['C']=[i for i, x in enumerate(RC) if x == 'C']
        location_dict['R']=[i for i, x in enumerate(RC) if x == 'R']
        C_div=[block_div_genome_1[i] for i in location_dict['C']]
        R_div=[block_div_genome_1[i] for i in location_dict['R']]
        if len(C_div)>0: 
            hmm_mu=sum(C_div)/len(C_div)
        if len(R_div)>0:
            hmm_mnb=sum(R_div)/len(R_div)
            merged_window_list=loc_to_pois(location_dict['R'], slide=50, winlen=window_size)
            ea=get_extracted_alignments(concat_S1S2[2],concat_S1S2[3], merged_window_list)
            cf=get_clonal_frame(concat_S1S2[2],concat_S1S2[3], merged_window_list,winlen=window_size)
#            if not os.path.exists(final_output_dir+'/alignment_fasta'):
#                os.mkdir(final_output_dir+'/alignment_fasta')
#            if not os.path.exists(final_output_dir+'/recomb_windows'):
#                os.mkdir(final_output_dir+'/recomb_windows')
#            if not os.path.exists(final_output_dir+'/clonal_frame'):
#                os.mkdir(final_output_dir+'/clonal_frame')
#            write_fasta(strain1=strain1,strain2=strain2,extracted_alignments=ea,merged_window_list=merged_window_list,final_output_dir=final_output_dir)
#            write_cf_fasta(strain1=strain1,strain2=strain2,extracted_alignments=cf,merged_window_list=merged_window_list,final_output_dir=final_output_dir, subfolder='/clonal_frame/')
#            write_csv(strain1,strain2,merged_window_list,final_output_dir=final_output_dir)

    block_div_genome_str = ','.join(str(v) for v in block_div_genome)
    edge = '\t'.join([strain1,
                     strain2,
                     str(init_div_1),
                     str(init_div),
                     str(alignment_size_1),
                     str(alignment_size),
                     str(g1size),
                     str(g2size),
                     str(observed_sum_sq_diff),
                     str(low_percentile),
                     str(high_percentile),
                     str(mu/window_size),
                     str(mnb/window_size),
                     str(100-fc),
                     str(a),
                     str(hmm_fr),
                     str(hmm_mu),
                     str(hmm_mnb)])
    

    return [edge, block_dg_1]

def logLL(log_par, data):

    mu,mnb,a,fc = log_par  #mnb is mean for nb,size=a
    p=a/(a+mnb)
    R1=poisson.pmf(data,mu)*fc
    R2=nbinom.pmf(data,a,p)*(100-fc)
    R1=np.nan_to_num(R1)
    R2=np.nan_to_num(R2)
    logL = sum(np.log(R1+R2+10**(-32)))
    return logL

def P_mle(data,block_div_genome,winlen):
    startmu=sum(block_div_genome)/len(block_div_genome)*winlen
    data=[d*winlen for d in data]
    bnds=((10**-5,winlen/10),(10**-5,winlen/3),(2,10),(0.1,99.9))
    res = opt.minimize(fun=lambda log_params, data: -logLL(log_params, data),
    x0=np.array([startmu,startmu*2,2,50]), args=(data,), method='L-BFGS-B',bounds=bnds,
    options={'ftol': 2.220446049250313e-07, 'gtol': 1e-04, 'maxls': 50})
                
    if res.success:
        if res.x[3]<=52 and res.x[3]>=48 and res.x[0]<=0.5:
  #          print (res)
            res=opt.minimize(fun=lambda log_params, data: -logLL(log_params, data),x0=np.array([startmu,startmu*2,2,50]), args=(data,), method='SLSQP',bounds=bnds)
  #          print (res)
        return res.x
    else:
        res=opt.minimize(fun=lambda log_params, data: -logLL(log_params, data),
    x0=np.array([startmu,startmu*2,2,50]), args=(data,), method='SLSQP',bounds=bnds)
        if res.success:
            return res.x
        else:
            return np.array([0,0,1,100])

def run_hmm(data,window_size,mle_output):
    #Set up model
    model = HiddenMarkovModel()
    #Set up parameters
    mu=mle_output[0]
    mnb=mle_output[1]
    a=mle_output[2]
    fc=mle_output[3]
    #Set up emission matrix
    maxSNP=max(np.array(data)*window_size)
    SNP=np.array(data)*window_size
    SNP=[str(round(x,0)) for x in SNP]
    SNP_x=numpy.arange(start=0,stop=maxSNP+1)
    dpois={}
    dnbinom={}
    p=a/(a+mnb)
    if mu<=10**-5+10**-6:
        mu=np.mean(np.array(data)*window_size)
    for x in SNP_x:
        dpois[str(x)]=poisson.pmf(x,mu)
        dnbinom[str(x)]=nbinom.pmf(x,a,p)
    
    dpois = DiscreteDistribution(dpois)
    dnbinom = DiscreteDistribution(dnbinom)
    spois = State(dpois, name="C")
    snbinom = State(dnbinom, name="R")
    model.add_states([spois, snbinom])
    #Set up transition matrix
    model.add_transition(model.start, spois, fc/100)
    model.add_transition(model.start, snbinom, (100-fc)/100)
    model.add_transition(spois, spois, 0.5)
    model.add_transition(snbinom, snbinom, 0.5)
    model.add_transition(spois, snbinom, 0.5)
    model.add_transition(snbinom, spois, 0.5)
    model.bake()
    
    #model.dense_transition_matrix()
    #Run viterbi learning
    model.fit([np.array(SNP)],algorithm='viterbi')
    model.bake()
    #model.dense_transition_matrix()
    RC=[]
    if model.viterbi(SNP)[1] is not None:
        RC=[state.name for i, state in model.viterbi(SNP)[1]]
    return RC

def loc_to_pois(R_list, slide, winlen):
    '''
    

    Parameters
    ----------
    R_list : list of window numbers classfied as R
    slide : step size for sliding windows
    winlen : length of window

    Returns: basepair positions for sequence extraction
    -------
    None.

    '''
    
    window_list=[]
    for w in R_list:
        start=(w*slide)-1
        end=start+winlen+1
    #    end=start+slide+1
        window_list.append([end,start])
    window_list=sorted(window_list,reverse=True)
    if len(window_list)>0:
        merged_window_list_r=merge_intervals(window_list)
        merged_window_list=[(x[1],x[0]-winlen) for x in merged_window_list_r]
        merged_window_list=sorted(merged_window_list)
    else:
        merged_window_list=[]
    return merged_window_list

def get_extracted_alignments(sequence_1, sequence_2, positions_to_remove):
    '''
    Helper method that extracts substrings when given a list of
    start and end positions to remove
    '''
    if positions_to_remove == []:
        return []
    else:
        merged = positions_to_remove
        starts, ends = zip(*sorted(merged))
    final_blocks=[]   
    for i, start_of_deleted_region in enumerate(starts):
        end_of_deleted_region = ends[i]
        subsequence_1 = sequence_1[start_of_deleted_region:end_of_deleted_region]
        subsequence_2 = sequence_2[start_of_deleted_region:end_of_deleted_region]
        subsequence_1 = subsequence_1.replace('0','')
        subsequence_1 = subsequence_1.replace('1','')
        subsequence_2 = subsequence_2.replace('0','')
        subsequence_2 = subsequence_2.replace('1','')
        final_blocks.append((subsequence_1, subsequence_2))
    return final_blocks

def get_clonal_frame(sequence_1, sequence_2, positions_to_remove,winlen):
    '''
    Helper method that extracts the clonal frame when given a list of
    start and end positions to remove
    To be conservative we remove extra 1 window length upstream and downstream of recombined fraction
    '''
    if positions_to_remove == []:
        return []
    else:
        merged = positions_to_remove
        starts, ends = zip(*sorted(merged))
        subsequence_1 = sequence_1[0:(starts[0]-winlen)]
        subsequence_2 = sequence_2[0:(starts[0]-winlen)]
        subsequence_1 = subsequence_1.replace('0','')
        subsequence_1 = subsequence_1.replace('1','')
        subsequence_2 = subsequence_2.replace('0','')
        subsequence_2 = subsequence_2.replace('1','')
        if starts[0]-winlen>0:
            final_blocks=[(subsequence_1, subsequence_2)]
        else:
            final_blocks=[]

    for i, end_of_deleted_region in enumerate(starts):
        if i == len(ends)-1:
            end_of_deleted_region = ends[i]+winlen
            start_of_deleted_region=len(sequence_1)+1
            subsequence_1 = sequence_1[end_of_deleted_region:start_of_deleted_region]
            subsequence_2 = sequence_2[end_of_deleted_region:start_of_deleted_region]
        else:
            start_of_deleted_region = starts[i+1]-winlen
            subsequence_1 = sequence_1[end_of_deleted_region:start_of_deleted_region]
            subsequence_2 = sequence_2[end_of_deleted_region:start_of_deleted_region]
        subsequence_1 = subsequence_1.replace('0','')
        subsequence_1 = subsequence_1.replace('1','')
        subsequence_2 = subsequence_2.replace('0','')
        subsequence_2 = subsequence_2.replace('1','')
        if start_of_deleted_region-end_of_deleted_region>0:
            final_blocks.append((subsequence_1, subsequence_2))
    return final_blocks

def write_fasta(strain1,strain2,extracted_alignments,merged_window_list,final_output_dir,
                subfolder='/alignment_fasta/'):
    outfile_1=final_output_dir + subfolder + strain1 +'-'+ strain2 +'-' + strain1 + '.fna'
    outfile_2=final_output_dir + subfolder + strain1 +'-'+ strain2 +'-' + strain2 + '.fna'
    ofile = open(outfile_1, "w")
    for i in range(len(extracted_alignments)):
        ofile.write(">" + strain1 + '-' + str(merged_window_list[i][0]+1) + '-' + str(merged_window_list[i][1]) + "\n" + extracted_alignments[i][0] + "\n")
    ofile.close()
    ofile = open(outfile_2, "w")
    for i in range(len(extracted_alignments)):
        ofile.write(">" + strain2 + '-' + str(merged_window_list[i][0]+1) + '-' + str(merged_window_list[i][1]) + "\n" + extracted_alignments[i][1] + "\n")
    ofile.close()

def write_cf_fasta(strain1,strain2,extracted_alignments,merged_window_list,final_output_dir,
                subfolder='/alignment_fasta/'):
    extracted_consensus=[]
    for i in range(len(extracted_alignments)):
        ea1=extracted_alignments[i][0]
        ea2=extracted_alignments[i][1]
        ea=''
        for s in range(len(ea1)):
            if ea1[s]==ea2[s]:
                ea=ea+ea1[s]
            else:
                ea=ea+random.choice([ea1[s],ea2[s]])
        extracted_consensus.append(ea)
    N500='N'*500
    extracted_consensus_seq=''
    for ea in extracted_consensus:
        extracted_consensus_seq=extracted_consensus_seq+ea+N500

    outfile=final_output_dir + subfolder + 'consensus_core_cf.fna'
    ofile = open(outfile, "w")
    ofile.write(">consensus_core_cf\n"+extracted_consensus_seq+"\n")
    ofile.close()

#    outfile_1=final_output_dir + subfolder + strain1 +'-'+ strain2 +'-' + strain1 + '.fna'
#    outfile_2=final_output_dir + subfolder + strain1 +'-'+ strain2 +'-' + strain2 + '.fna'
#    ofile = open(outfile_1, "w")
#    for i in range(len(extracted_alignments)):
#        ofile.write(">" + strain1 + '-' + str(i) + "\n" + extracted_alignments[i][0] + "\n")
#    ofile.close()
#    ofile = open(outfile_2, "w")
#    for i in range(len(extracted_alignments)):
#        ofile.write(">" + strain2 + '-' + str(i) + "\n" + extracted_alignments[i][1] + "\n")
#    ofile.close()


def write_csv(strain1,strain2,merged_window_list,final_output_dir):
    df = pd.DataFrame.from_records(merged_window_list, columns =['start', 'end'])
    df['start']=df['start']+1  #this is for changing into a 'normal' expression instead of python index
    
    df.to_csv(final_output_dir + '/recomb_windows/' + strain1 + '-' + strain2 + '-recombwin.csv')

def get_cumulative_window_spectrum(idw, gs):
    '''
    Gets the X and Y coordinates of the identical window spectrum
    i.e., the fraction of the genome belonging to identical windows
    above a certain size
    '''

    obs_frac_counts = np.zeros(gs)
    norm = np.sum(idw)
    windows = Counter(idw)
    for wsize, count in windows.items():
        obs_frac_counts[wsize] = count * wsize * 1.0 / norm
    return 1.0 - np.cumsum(obs_frac_counts)


def get_concatenated_alignment(alignment):
    '''
    This creates a list of tuples that constitute a concatenated alignment.
    Every entry in the list is a tuple that corresponds to an alignment block.
    '''

    with open(alignment, 'r') as infile:
        '''
        Parser assumes a maf format where every alignment block begins with a
        statement of how many sequences are in that block, indicated by
        "mult=." Also assumes that the order of sequences in each block is
        the same.
        '''
        seqs = []
        total_len = 0
        for lines in infile:
            if 'mult=2' in lines:
                seq_line_1 = next(infile)
                block_1 = seq_line_1.split()[-1].strip()
                total_len += len(block_1)
                seq_line_2 = next(infile)
                block_2 = seq_line_2.split()[-1].replace('\n', '')
                seqs.append((block_1, block_2))
    return seqs, total_len


def id_var_window_counts(sequence_1, sequence_2):
    '''
    This method takes two aligned sequences (strings) and returns the
    lengths of all identical windows between those sequences.
    '''
    if sequence_1 == sequence_2:
        id_seqs = [len(sequence_1)]
    else:
        a1 = np.array(list(sequence_1))
        a2 = np.array(list(sequence_2))
        mutated_positions = np.where(a1 != a2)[0]
        id_seqs = -1 + np.ediff1d(mutated_positions,
                                  to_begin=mutated_positions[0] + 1,
                                  to_end=len(sequence_1) - mutated_positions[-1])
    return id_seqs


def naive_div_count(sequence_1, sequence_2):
    '''
    Given two aligned strings, returns the number of differences between them.
    '''
    if sequence_1 == sequence_2:
        return 0
    a1 = np.array(list(sequence_1))
    a2 = np.array(list(sequence_2))
    return len(np.where(a1 != a2)[0])

def block_div_calc(sequence_1, sequence_2, winlen, slide=50):
    '''
    Calculate diversity in blocks across a genome
    :param sequence_1:
    :param sequence_2:
    :param winlen:
    :return:
    '''
    winlen=int(winlen)
    begin = 0
    div_list = []
    for end in range(winlen, (len(sequence_1)-winlen), slide):
        d = naive_div_count(sequence_1[begin:(begin+winlen)], sequence_2[begin:(begin+winlen)])
        div_list.append(d/winlen)
        begin = end

    if begin < len(sequence_1):
        d = naive_div_count(sequence_1[begin:], sequence_2[begin:])
        div_list.append(d/len(sequence_1[begin:]))
    return div_list

def filter_block(sequence_1, sequence_2):
    removal_positions = filter_string(sequence_1)
    removal_positions += filter_string(sequence_2)
    return [block for block in get_filtered_subblocks(sequence_1, sequence_2, removal_positions) if block[0] != '']


def filter_string(S):
    groups = groupby(S)
    result = [(label, sum(1 for _ in group)) for label, group in groups]
    begin = 0
    filter_intervals = []
    for base, count in result:
        end = begin + count
        if base == '-' and count >= 2:
            filter_intervals.append((end, begin))
        begin += count
    return(filter_intervals)


def filter_block_by_divergence(sequence_1, sequence_2, init_div, winlen, winnumber,binom_cutoff,slide):
    '''
    Filters two sequences from an alignment block to remove regions
    that are significantly more diverged than expected
    '''
    if sequence_1 == sequence_2:
        return [(sequence_1, sequence_2)]

    else:
        removal_positions = []
        begin = 0
        end=begin+winlen
        
        while end< len(sequence_1):        
            d = naive_div_count(sequence_1[begin:end], sequence_2[begin:end])
            if d >= binom_cutoff:
                removal_positions.append((end, begin))
               
            begin = begin+slide
            end = begin+winlen

        if begin < len(sequence_1):
            d = naive_div_count(sequence_1[begin:], sequence_2[begin:])
            binom_cutoff_2 = binom.ppf(q=1-0.05/winnumber,n=(len(sequence_1) - begin),p=init_div,loc=0)
            if d >= binom_cutoff_2:
                removal_positions.append((len(sequence_1), begin))
        return [block for block in get_filtered_subblocks(sequence_1, sequence_2, removal_positions) if block[0] != '']


def get_filtered_subblocks(sequence_1, sequence_2, positions_to_remove):
    '''
    Helper method that splits a string into substrings when given a list of
    start and end positions to remove
    '''
    if positions_to_remove == []:
        return [(sequence_1, sequence_2)]
    else:
        final_blocks = []
        initial_start = 0
        if len(positions_to_remove) > 1:
            merged = merge_intervals(sorted(positions_to_remove, reverse=True))
        else:
            merged = positions_to_remove
        ends, starts = zip(*sorted(merged, key=lambda x: x[1]))
        for i, start_of_deleted_region in enumerate(starts):
            end_of_deleted_region = ends[i]
            subsequence_1 = sequence_1[initial_start:start_of_deleted_region]
            subsequence_2 = sequence_2[initial_start:start_of_deleted_region]
            initial_start = end_of_deleted_region
            final_blocks.append((subsequence_1, subsequence_2))
        final_blocks.append((sequence_1[initial_start:], sequence_2[initial_start:]))
        return final_blocks


def single_param_null_model(w, div):
    '''
    The simple single parameter null model that describes
    the window spectrum under an assumption of only mutation
    and no transfer
    '''
    return np.exp(-div * w) * (div * w + 1)


def merge_intervals(intervals):
    all_intervals = []
    for j, interval in enumerate(intervals):
        end, start = interval
#        print (interval)
        if j == 0:
            current_interval = interval
        else:
            if intervals[j-1][1] <= end:
                current_interval = (current_interval[0], start)
            else:
                all_intervals.append(current_interval)
                current_interval = interval
#        print (current_interval)
    if len(all_intervals) > 0:
        if all_intervals[-1] != current_interval:
            all_intervals.append(current_interval)
    else:
        all_intervals.append(current_interval)
    return all_intervals


