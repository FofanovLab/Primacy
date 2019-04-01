import pandas as pd
import numpy as np
from multiprocessing import Pool
from Bio.Seq import Seq


from primacy.utils import (
    get_json_obj, write_json_obj,
    convert_dict_to_dataframe,
    get_primer_dataframe)
from primacy.primer import get_existing_primer_stats, specificity
from primacy.sequence import Sequence

def get_random_primer(flank):
    return np.random.choice(flank)

def get_random_additional_values():
    return {
        "d_tm": float(np.random.normal(0, 5)),
        "d_avg_tm": float(np.random.normal(0, 10)),
        "d_avg_gc": float(np.random.normal(0,5)),
        "avg_cross_hybrid_score": float(np.random.uniform(1,50)),
        "avg_cross_hybrid_percent": float(np.random.uniform(0, 100)),
        "avg_cross_hybrid_run": int(np.random.choice(range(0, 10))),
        "max_cross_hybrid_score": float(np.random.uniform(1,50)),
        "max_cross_hybrid_percent": float(np.random.uniform(0, 100)),
        "max_cross_hybrid_run": int(np.random.choice(range(0,10)))}
        
def get_random_additional_pos_values():
    return {"amplicon_size": int(np.random.choice(range(100, 300))),
        "targ_dist_from_right_primer": int(np.random.choice(range(0,50))),
        "targ_dist_from_left_primer": int(np.random.choice(range(0, 50)))}

def get_optimized_primer_sets(
    file_paths, workdir, opt_params,
    pc_params, threads):
    optimized = {}
    for file_path in file_paths:
        primer_obj = get_json_obj(file_path)
        for seq_id, data in primer_obj.items():
            optimized[seq_id] = {"forward": {}, "reverse": {}}
            if opt_params['include'][seq_id]["forward"]:
                forward_id = get_random_primer(
                    opt_params['include'][seq_id]["forward"])
                optimized[seq_id]['forward'][forward_id] = \
                    {**data['forward'][forward_id],
                    **get_random_additional_pos_values(),
                    **get_random_additional_values()}
            if opt_params['include'][seq_id]["reverse"]:
                reverse_id = get_random_primer(
                    opt_params['include'][seq_id]["reverse"])
                optimized[seq_id]['reverse'][reverse_id] = \
                    {**data['reverse'][reverse_id],
                    **get_random_additional_pos_values(),
                    **get_random_additional_values()}

    opt_params['result'] = optimized

            
    # for primer_id, seq in params['background'].items():
    df = pd.DataFrame.from_dict(
        opt_params['params']['background'], orient='index')
    df.columns = ["sequence"]
    df['sequence'] = df.sequence.apply(
        lambda x: Seq(x, Sequence.alphabet).upper() )
    df = specificity(
        df,
        pc_params['background_seq'],
        workdir, threads)
    background_primers = get_existing_primer_stats(
            df, pc_params["pcr_salts"]).to_dict(orient='index')
    opt_params['background_result'] = {
        **background_primers,
        **get_random_additional_values()}
    opt_params['set_score'] = 0.8
    return opt_params
    



