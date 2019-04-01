import logging
import subprocess
import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp
from itertools import groupby, compress, product
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

AMB = IUPAC.IUPACData.ambiguous_dna_values
HYBRID_SCORES = {'A': 2, 'T': 2, 'C':4, 'G': 4}

def expand_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = IUPAC.IUPACData.ambiguous_dna_values
    return tuple(map("".join, product(*map(d.get, seq))))

def ambiguous_count(seq):
    """count the number of possible sequences from ambiguous DNA input"""
    d = IUPAC.IUPACData.ambiguous_dna_values
    return np.prod([len(d[base]) for base in seq])


def gc(primer):
    """Returns the percent of GC content in the primer. This will
    count any ambiguous base that specifies either a G or G but 
    will ignore any other ambiguous base.
    """
    return GC(primer)


def gc_clamp(primer):
    """
    Returns the number (int) of G or C's in the last 5 bases (3' end)
    of the primer. This includes ambiguous bases that specify 
    either a G or C but ignores any other ambiguous bases (e.g. N)
    """
    # get length in (unlikely) case where primer is shorter than 5
    return int(GC(primer[-5:]) * len(primer[-5:]) / 100)

def tm(tm_params):
    def tm_func(primers):
        """ Calculates the nearest neighbor melting temperature using
        the user specified pcr salt parameters. When there are multiple
        sequences due to ambiguous bases, an average Tm is returned.
        """
        return np.mean([MeltingTemp.Tm_NN(
            primer, **tm_params) for primer in primers])
    return tm_func


def homopolymers(primer_combos):
    """Calculates the number of bases that are in homopolymer
    runs greater than length 3 and determines the longest run.
    Returns a dictionary containing the percent of bases that are
    in a homopolymer stretch and the length of the longest run.
    For primers with ambiguous bases, it returns the worst case for
    each which may belong to different primer sequences.
    """
    base_counts = []
    longest = []
    for primer in primer_combos:
        homopolymers = [list(hp) for run, hp in groupby(str(primer))]
        bases = []
        for hp in homopolymers:
            if len(hp) > 3:
                bases.append(len(hp))
        base_counts.append(sum(bases))
        longest.append(max(bases, default=0))
    return {
        'percent': (max(base_counts) / len(primer)) * 100,
        'run': max(longest)}


def dimerization(primer_1, primer_2):
    """
    Returns a dictionary with the percent of dimerization and
    the length of the longest run (weighted by C-G vs. A-T pairings),
    the median length of complementary runs, and an overhang score. 
    The values are each the maximum (worst-case) possiblities
    among all alignments of the primers. Provides 
    worst-case for primers with ambiguous bases.
    primer_1 (str, Seq): First primer in 5' to 3' direction
    primer_2 (str, Seq): Second primer reverse_complemented
                         in 3' to 5' direction.
    """
    p1 = primer_1 + ("*" * (len(primer_2) - 1))
    p2 = "*" * (len(primer_1) - 1) + primer_2[::-1]
    max_percent = []
    max_run = []
    max_median = []
    max_overhang_run = []
    while len(p1):
        filter_ = [0 if "*" in primers else 1 for primers in zip(p1, p2)]
        percent, run, median, complement = dimerization_worker(
            compress(p1, filter_),
            compress(p2, filter_))
        max_run.append(run)
        max_percent.append(percent)
        max_median.append(median)
        if "*" in p1:
            # "*" indicates primer is in 5' overhang state
            max_overhang_run.append(sum(complement))
        p1 = p1[:-1]
        p2 = p2[1:]    
    return {
        'percent': (max(max_percent)/(max([len(primer_1), len(primer_2)]) * 4)) * 100,
        'run': max(max_run),
        'median': max(max_median),
        'hairpin': max(max_overhang_run)}
  

def dimerization_worker(primer_1, primer_2):
    """ Returns the total number of complementary bases and the longest
    run of complementary bases (weighted by HYBRID_SCORES), the median
    length of all runs and the array of complementary bases. 
    """
    p1 = [set(AMB[i]) for i in primer_1]
    p2 = [set(AMB[i]) for i in primer_2]
    complementary = [
        max(
            [HYBRID_SCORES[base] for base in s1.intersection(s2)],
            default=0) for s1, s2 in zip(p1, p2)]
    total = sum(complementary)
    complementary_runs = [list(comp) for run, comp in groupby(
        complementary, key=lambda x: x > 0)]
    max_run = sum(max(complementary_runs, key=sum))
    run_lens = [sum(c) for c in complementary_runs if sum(c) > 0]
    if run_lens:
        run_lens.append(0)
    median = np.median(run_lens) if run_lens else 0
    return (total, max_run, median, complementary)

def get_existing_primer_stats(df, tm_params):
    df['ambiguous_seqs'] = df.sequence.apply(
        expand_ambiguous_dna
    )
    df['degenerate'] = df.ambiguous_seqs.apply(
        lambda x: len(x)
    )
    df['gc'] = df.sequence.apply(gc)
    df['gc_clamp'] = df.sequence.apply(gc_clamp)
    tm_calc = tm(tm_params)
    df['tm'] = df.ambiguous_seqs.apply(tm_calc)
    df['homopolymers'] = df.ambiguous_seqs.apply(homopolymers)
    df['reverse_complement'] = df.sequence.apply(
        lambda x: x.reverse_complement())
    df['dimerization'] = df[['sequence', 'reverse_complement']].apply(
        lambda row: dimerization(
            row.sequence, row.reverse_complement), axis=1)
    df["seq"] = df.sequence.astype(str)
    df.drop(
        ["ambiguous_seqs", "reverse_complement", "sequence"],
        axis=1, inplace=True)
    return df


def get_primer_stats(df, tm_params):
    df['primer_id'] = df[[
        "seq_id", "start", "length", "flank"]].apply(
            lambda row: "{0}_{1}_{2}_{3}".format(
                row.seq_id, row.start, row.length, row.flank), axis=1)
    df["flank"] = df.flank.apply(lambda x: 0 if x=="F" else 1)

    return get_existing_primer_stats(df, tm_params)



def specificity(primer_df, background_paths, outpath, threads):
    primer_df['specificity'] = [
        {'mm0': 0, 'mm1': 0, 'mm2': 0, 'mm3': 0} for _ in primer_df.index]
    return primer_df


# def specificity_setup(primer_df, background_paths, outpath, threads):
#     primer_df['specificity'] = {'mm0': 0, 'mm1': 0, 'mm2': 0, 'mm3': 0}
#     background_fasta = os.path.join(outpath, "background.fasta")
#     primer_query_fasta = os.path.join(outpath, "primer_queries.fasta")
#     database_path = os.path.join(outpath, "blast_db/blast_db")
#     blast_out = os.path.join(outpath, "blast_db/blast.xml")
#     n_seqs = get_background_sequences(
#         background_paths, background_fasta)
#     get_primer_fasta(
#         primer_df.sequence, primer_df.primer_ids, primer_query_fasta)
#     try:
#         make_blast_db(background_fasta, database_path)
#     except subprocess.CalledProcessError:
#         # just return 0 values
#         return primer_df
#     run_blast(
#         database_path, primer_query_fasta, blast_out, threads,
#         **{"word_size": 7, "max_target_seqs": n_seqs}
#     )
    

# def run_blast(
#     db_path, primer_query_fasta, outpath_xml, threads, **qwargs):
#     blastx_cline = NcbiblastnCommandline(
#         query=primer_query_fasta, db=db_path,
#         outfmt=5, out=outpath_xml,
#         num_threads=threads)
#     blastx_cline()
#     result_handle = open(outpath_xml)
#     blast_records = NCBIXML.parse(result_handle)
#     for blast_record in blast_records:
#         pass




def get_background_sequences(background_paths, background_fasta):
    """
    Combines each of the background fasta files into the same fasta
    file: background fasta.
    background_paths (list): list of paths to fasta files
    background_fasta (str): output file name
    Returns the total number of sequences
    """
    total_sequences = 0
    for file_path in background_paths:
        with open(file_path, 'r') as fasta_handle:
            for record in SeqIO.parse(fasta_handle, 'fasta'):
                total_sequences += 1
                SeqIO.write(record, background_fasta, 'fasta')
    return total_sequences



def get_primer_fasta(primers, primer_ids, fasta_output):
    """
    Writes each primer sequence to the fasta_output file in 
    fasta format with the primer_id and the header
    """
    with open(fasta_output, 'w') as fasta:
        for primer_id, primer in zip(primer_ids, primers):
            SeqIO.write(
                SeqRecord(primer, id=primer_id), fasta, 'fasta')

    

def make_blast_db(subject_fp, db_path):
    """
    Makes a subprocess call to build blast database. 
    subject_fp (str): path to background sequence fastas file
    db_path (str): path to database
    Raises CalledProcessError if subprocess fails.
    """
    command = ["makeblastdb", "-in", subject_fp, "-out", db_path,
               "-dbtype", "nucl"]

    try:
        subprocess.call(command)
    except subprocess.CalledProcessError:
        logger.warn(
            "Unable to generate blast database for specificity "
              "alignment calculation")
        raise
