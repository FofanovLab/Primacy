import json
import pandas as pd
# from dask import dataframe as dd
# from dask.multiprocessing import get
 

def read_config(config):
    with open(config, 'r') as json_in:
        return json.loads(json_in.read())

def get_json_obj(json_path):
    return read_config(json_path)

def write_json_obj(json_obj, json_path):
    with open(json_path, 'w') as json_out:
        json.dump(json_obj, json_out)

def get_primer_obj_paths(config_obj):
    return [
        seq_obj['outfile'] for seq_obj
        in config_obj['sequences'].values()]


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


# def parallel_apply(df, func, column_map, threads, meta):
#     return dd.from_pandas(
#         df, npartitions=threads).\
#             map_partitions(
#             lambda df : df.apply(
#          lambda x : func(**{
#              k: x[v] for k, v in column_map.items()}), axis=1),
#              meta=meta).\
#             compute(scheduler='threads')



