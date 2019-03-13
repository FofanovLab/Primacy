import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
from primacy.utils import get_json_obj, write_json_obj

def min_max_scale_clipped(values, v_min, v_max):
    scaled = np.divide(np.subtract(values, v_min), np.subtract(v_max, v_min))
    scaled[scaled > 1] = 1
    return scaled

def get_distance(values, range_min, range_max):
    dist = np.zeros(len(values))
    above = np.where(values > range_max)
    below = np.where(values < range_min)
    dist[above] = np.abs(np.subtract(values[above], range_max))
    dist[below] = np.abs(np.subtract(range_min, values[below]))
    return dist


def convert_dict_to_dataframe(primer_obj):
    """
    Combines forward and reverse primers from primer_obj into a data frame
    """
    return pd.DataFrame.from_dict(
        {**primer_obj['forward'], **primer_obj['reverse']},
        orient="index")


def get_primer_dataframe(df):
    """
    Return a dataframe with expanded dictionary
    columns (e.g. homopolymer, specificity, and dimerization)
    into individual columns.
    """
    # expand dictionary columns into separate columns
    return pd.concat(
            [df.drop(
                ['specificity', 'homopolymers', 'dimerization'], axis=1),
                df['specificity'].apply(pd.Series),
                df['homopolymers'].apply(pd.Series).rename(
                    columns={
                        c: "homo_{}".format(c) for c in ["percent", "run"]}),
                df['dimerization'].apply(pd.Series).rename(
                    columns={c: "dimer_{}".format(c) for c in [
                        "percent", "run", "median", "hairpin"]})
                ], axis=1)


def convert_to_distance(primer_df, tm_opt, gc_opt, gc_clamp_opt=2):
    """
    Convert tm, gc%, and gc_clamp to an absolute distance
    (tm_dist, gc_dist, gc_clamp_dist) 
    away from optimum range. This makes it so that all features will need
    to be minimized. 
    """
    primer_df['tm_dist'] = get_distance(
        primer_df.tm.values, tm_opt, tm_opt)
    primer_df['gc_dist'] = get_distance(
        primer_df.gc.values, gc_opt['min'], gc_opt['max'])
    primer_df['gc_clamp_dist'] = get_distance(
        primer_df.gc_clamp.values, gc_clamp_opt, gc_clamp_opt)
    # primer_df.drop(['tm', 'gc', 'gc_clamp'], axis=1, inplace=True)
    return primer_df


def scale_weights(weights, index):
    # Some weights are broken down into muliple components 
    # e.g. mm0, mm1, mm2, mm3 for specificity. So the weight
    # adj determines the split
    weight_adj = {
        "degenerate": {"degenerate": 1},
        "gc": {"gc_dist": 0.7, "gc_clamp_dist": 0.3},
        "tm": {"tm_dist": 1},
        "specificity": {"mm0": 0.75, "mm1": 0.2, "mm2": 0.04, "mm3": 0.01},
        "homopolymer": {"homo_percent": 0.5, "homo_run": 0.5},
        "dimer": {"dimer_run": 0.35, "dimer_percent": 0.35,
            "dimer_hairpin": 0.20, "dimer_median": 0.10}}
    adj_weights = pd.Series(np.zeros(len(index)), index=index)
    denominator = np.sum(np.array(list(weights.values()), dtype=float))
    for group, weight in weights.items():
        for sub_group, adjustment in weight_adj[group].items():
            adj_weights[sub_group] = (weight / denominator) *  adjustment
    return adj_weights.values


def get_max_for_scale(index):
    max_dic = {"degenerate": 1024,
                "gc_dist": 30, "gc_clamp_dist": 3,
                "tm_dist": 30, "mm0": 3, "mm1": 5, "mm2": 10, "mm3": 20,
                "homo_percent": 50, "homo_run": 6, "dimer_percent": 50,
                "dimer_run": 28, "dimer_median": 16, "dimer_hairpin": 16}
    return pd.Series(
        [max_dic[i] for i in index], index=index).values.reshape(
            len(index), 1).T

def get_min_for_scale(index):
    scale_min = pd.Series(np.zeros(len(index)), index=index)
    scale_min['degenerate'] = 1
    return scale_min.values.reshape(len(index), 1).T

def apply_weights(scaled_matrix, weights):
    return scaled_matrix.dot(weights.T)

def get_scores(primer_df, weights):
    index = [
    "degenerate", "gc_dist", "gc_clamp_dist", "tm_dist",
    "mm0", "mm1", "mm2", "mm3", "homo_percent", "homo_run",
    "dimer_percent", "dimer_run", "dimer_median", "dimer_hairpin"]
    scale_min = get_min_for_scale(index)
    scale_max = get_max_for_scale(index)
    adj_weights = scale_weights(weights, index)
    primer_df['score'] = apply_weights(min_max_scale_clipped(
        primer_df[index].values, scale_min, scale_max), adj_weights)
    # add a percent rank where smaller percentile is better score
    primer_df['rank'] = primer_df.groupby('flank')['score'].rank(
        ascending=True, pct=True)
    return primer_df


def add_scores_to_primer_obj_helper(row, primer_obj, flank):
    primer_obj[flank][row.name]['rank'] = row['rank']
    primer_obj[flank][row.name]['score'] = row['score']



def add_scores_to_primer_obj(primer_obj, primer_df):
    primer_df[primer_df.flank == 0].apply(
        lambda row: add_scores_to_primer_obj_helper(
            row, primer_obj, "forward"), axis=1 )    
    primer_df[primer_df.flank == 1 ].apply(
        lambda row: add_scores_to_primer_obj_helper(
            row, primer_obj, "reverse"), axis=1 )
    return primer_obj


def get_primer_score(primer_obj, tm_opt, gc_opt, weights):
    """
    Returns a primer_obj with scores and ranks added.
    primer_obj (dict): {"forward": {"primer_id": {...}, ...},
                        "reverse": {"primer_id": {...}, ...}}
    tm_opt (int): optimal melting temp for PCR
    gc_opt (dict): optimal range for GC content {"min": 30, "max": 60}
    weights (dict): relative importance for each feature group
                    {<feature name> : 1, ... }
    """
    # make dataframe from primer_obj, 
    # convert values (tm, gc, etc) to absolute distances away from opt so that
    # all features should be minimized. scale values and apply
    # weights to get a score and rank column, then add scores and ranks to 
    # primer object
    
    return add_scores_to_primer_obj(
        primer_obj,
        get_scores(
            convert_to_distance(
                get_primer_dataframe(convert_dict_to_dataframe(primer_obj)),
                tm_opt, gc_opt), weights))


def get_primer_scores_helper(file_path, params):
    primer_obj = get_json_obj(file_path)
    seq_id = list(primer_obj.keys())[0]
    return {seq_id: get_primer_score(
        primer_obj[seq_id], **params)}
    

def get_primer_scores_mp(file_paths, params, threads):
    p = Pool(threads)
    func = partial(get_primer_scores_helper, params=params)
    for file_path, primer_dic in zip(file_paths, p.imap(func, file_paths)):
        write_json_obj(primer_dic, file_path)
    p.close()
    p.join()






