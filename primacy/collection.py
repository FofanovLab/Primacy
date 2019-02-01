import json
import os
import pandas as pd
import logging
from multiprocessing import Pool
from functools import partial
from primer import get_primer_stats, specificity

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def remove_duplicate_primers(df):
    """ 
    Remove duplicate primers so that they do not 
    interfere with other amplicons
    """
    before = df.shape[0]
    df.drop_duplicates(subset="sequence", inplace=True)
    after = df.shape[0]
    logger.info("{} duplicate primers removed.".format(
        before-after))

def get_primer_df(seq_obj_list):
    """ 
    Get list of possible primers for sequence objects
    and add them to a dataframe with columns: seq_id, pos,
    length, flank, and sequence.
    """
    all_primers = []
    for seq_obj in seq_obj_list:
        all_primers += seq_obj.primers
    return pd.DataFrame(
        all_primers,
        columns=["seq_id", "start",
        "length", "flank", "sequence"])



def get_primer_data(seq_id_df, tm_params):
    """
    Calculate primer stats by seq_id and flank and 
    output a dictionary of values
    """
    seq_id, df = seq_id_df
    df = get_primer_stats(df, tm_params)
    primer_dict = {'forward': {}, 'reverse': {}}
    for flank, flank_df in df.groupby(['flank'], sort=False):
        flank_key = 'forward' if flank==0 else 'reverse'
        flank_df.index = flank_df.primer_id
        flank_df.drop(["primer_id"], inplace=True, axis=1)
        primer_dict[flank_key] = flank_df.to_dict(orient='index')
    return {seq_id: primer_dict}

def get_primer_collection(
        seq_obj_dict, tm_params, background_paths, outpath, threads):
    """
    Gathers upstream and downstream primers for each sequence object
    and calculates statistics for each primer. An output json is
    written for each sequence object named seq_<SEQ_ID>.json.
    seq_obj_dic: (dict) a mapping of seq_id to sequence objects
    tm_params: (dict) a dictionary of pcr salt parameters used to
            modify melting temperature calculation.
            ({"Na": 50, "K": 0, "Tris": 0, "Mg": 0, "dNTPs": 0})
    background_paths: (list) as list of paths to background sequences
            used to determine primer specificity.
    outpath: (str) path to write output files
    threads: (int) number of threads for multiprocessing
    """
    all_primers = get_primer_df(seq_obj_dict.values())
    remove_duplicate_primers(all_primers)
    all_primers = specificity(
        all_primers, background_paths, outpath, threads)
    p = Pool(threads)
    func = partial(get_primer_data, tm_params=tm_params)
    for primer_dic in p.imap(
        func, all_primers.groupby(["seq_id"], sort=False)):
        file_name = os.path.join(
                outpath,
                "seq_{}.json".format(list(primer_dic.keys())[0]))
        with open(file_name, 'w') as outfile:
            json.dump(primer_dic, outfile)
        



if __name__ == '__main__':
    pass