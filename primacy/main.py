import click
import logging
import json
import os

from primacy.sequence import Sequence
from primacy.collection import get_primer_collection
from primacy.primer_score import get_primer_scores_mp
from primacy.utils import read_config, get_json_obj, write_json_obj
# from primacy.score_existing_primers import score_existing_primers

@click.group()
def cli():
    pass

# @cli.command()
# @click.argument('project_name', type=str)
# @click.argument('fasta_path', type=click.File())
# @click.argument('outpath', type=click.Path(exists=False), default=os.getcwd())
# @click.argument('background_paths', nargs=-1)
# @click.option('--tm_opt', type=int, default=55)
# @click.option('--gc_min', type=int, default=30)
# @click.option('--gc_max', type=int, default=60)
# @click.option('--primer_weights', nargs=6,  default=[1,1,1,1,1,1])
# @click.option('--set_weights', nargs=3, default=[1,1,1])
# @click.option('--na', type=float, default=50)
# @click.option('--k', type=float, default=0)
# @click.option('--tris', type=float, default=0)
# @click.option('--mg', type=float, default=0)
# @click.option('--dntps', type=float, default=0)
# @click.option('--threads', type=int, default=2)
# def score_existing_primer_set(project_name, fasta_path, outpath,
#     background_paths, primer_weights, set_weights,
#     tm_opt=55, gc_min=30, gc_max=60,
#     na=50, k=0, tris=0, mg=0, dntps=0, threads=2):
#     tm_params = {'Na': na, "K": k, "Tris": tris, 'Mg': mg, "dNTPs": dntps}
#     gc_opt = {'min': gc_min, 'max': gc_max}
#     score_existing_primers(
#         project_name, fasta_path, tm_params, tm_opt, gc_opt, primer_weights,
#         set_weights, background_paths, outpath, threads)



@click.group()
def gui():
    pass

@gui.command()
@click.argument('config', type=click.Path(exists=True))
@click.argument('outpath', type=click.Path(exists=False), default=os.getcwd())
@click.option('--threads', '-t', default=2, type=click.IntRange(min=1))
def primer_collection(config, threads, outpath=os.getcwd()):
    config_obj = read_config(config)
    background_paths = \
        config_obj['primer_collection']['params']['background_seq']
    tm_params = config_obj['primer_collection']['params']['pcr_salts']
    seq_obj_dict = {seq_id: Sequence(
        seq_id=seq_id,
        primer_len_min=seq_info['primer_len_range']['min'],
        primer_len_max=seq_info['primer_len_range']['max'],
        **seq_info)
        for seq_id, seq_info in config_obj['sequences'].items()}
    # Run primer collection and return dictionary of output files
    output = get_primer_collection(
        seq_obj_dict, tm_params, background_paths, outpath, threads)
    # modify config json file with path to primer_collection for each sequence
    for seq_id, outpath in output.items():
        config_obj['sequences'][seq_id]['outfile'] = outpath
    write_json_obj(config_obj, config)


@gui.command()
@click.argument('config', type=click.Path(exists=True))
@click.option('--threads', '-t', default=2, type=click.IntRange(min=1))
def primer_score(config, threads):
    config_obj = read_config(config)
    primer_obj_paths = [
        seq_obj['outfile'] for 
            seq_obj in config_obj['sequences'].values()]
    get_primer_scores_mp(
        primer_obj_paths, config_obj['primer_scores']['params'], threads)


if __name__ == '__main__':
    gui()