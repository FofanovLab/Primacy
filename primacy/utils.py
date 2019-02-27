import json
from dask import dataframe as dd
from dask.multiprocessing import get


def read_config(config):
    with open(config, 'r') as json_in:
        return json.loads(json_in.read())

def get_json_obj(json_path):
    return read_config(json_path)

def write_json_obj(json_obj, json_path):
    with open(json_path, 'w') as json_out:
        json.dump(json_obj, json_out)

def parallel_apply(df, func, column_map, threads, meta):
    return dd.from_pandas(
        df, npartitions=threads).\
            map_partitions(
            lambda df : df.apply(
         lambda x : func(**{
             k: x[v] for k, v in column_map.items()}), axis=1),
             meta=meta).\
            compute(scheduler='threads')



