import subprocess
import os
import signal
import glob
import re
import sqlite3                                                                                                                                         
from sqlite3 import Error
import collections 

CommandResults = collections.namedtuple('CommandResults',['file_name',
                                                          'std_out',
                                                          'std_err'])

def run_command(command):
    return subprocess.Popen(command,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

def run_conditional(command, run_msg, dest_file, run_truth, write_truth):
    if not run_truth:
        print(run_msg)
        stds = run_command(command)
        if write_truth:
            write_file(stds[0].decode('utf-8'), dest_file)
    else:
        print('exists, not doing action: ' + run_msg)
        stds = (b'',b'')
    return CommandResults(dest_file,
                          stds[0].decode('utf-8'),
                          stds[1].decode('utf-8'))

def concatenate_files(file_list, dest_file):
    file_list_expanded = list(map(lambda x: x.file_name, file_list))
    command = ['./scripts/general/concatenate_files.sh'] + file_list_expanded
    run_msg='concatenating'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)


def convert_headers(source_file, dest_file, label):
    command = ['./scripts/canon/convert_headers.sh',
               source_file,
               label]
    run_msg='converting headers'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)
 
def blastp_multi(query ,source_file,dest_file):
    command = ['./scripts/nrdb/blast_nr.sh',
               source_file, query]
    run_msg='blasting multiple'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def make_uba_fasta(directory, search_string, dest_file):
    command = ['./scripts/uba/make_uba_fasta_db_direct.sh',
               directory, search_string]
    run_msg='making uba_fasta'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def make_index_map(source_file, dest_file):
    command = ['./scripts/canon/make_index_map.sh',
               source_file]
    run_msg='making index map'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def blast_map_to_fasta(source_file, label):
    command = ['./scripts/general/blast_to_fasta.sh',
               source_file.file_name, label]
    dest_file=source_file.file_name +'.fasta'
    run_msg='running blast map to fasta'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def make_blastp_db(source_file, out_database_name):
    command = ['makeblastdb', '-in', source_file, '-input_type', 'fasta',
               '-dbtype', 'prot', '-out', out_database_name]
    exists = (os.path.isfile(out_database_name+".pin") or
                   os.path.isfile(out_database_name+".phr") or
                   os.path.isfile(out_database_name+".psq"))
    run_msg='making blast db'
    return run_conditional(command, run_msg, out_database_name, exists, False)
    
def blastp(query, source_file, dest_file_name, evalue, coverage):
    command = ['blastp', '-query', query ,'-db', source_file,
               '-outfmt', '6 sseqid sseq evalue bitscore qcovsq',
               '-evalue', evalue,'-qcov_hsp_perc', coverage]
    run_msg='running blast map to fasta'
    exists=os.path.isfile(dest_file_name)
    return run_conditional(command, run_msg, dest_file_name, exists, True)

        
def get_max_info_blast_rep(blast_result_file):
    command = ['./scripts/general/max_blast_results_with_index.sh', blast_result_file.file_name]
    dest_file=blast_result_file.file_name+'.map'
    run_msg='running blast map to fasta'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)
    
def get_rgi_card_results(directory, search_string):
    file_name = 'tmp/results.ubaindex'
    command = ['./scripts/make_uba_index.sh', directory, search_string]
    run_msg='rcollecting card results'
    exists=os.path.isfile(file_name)
    return run_conditional(command, run_msg, file_name, exists, True)

def make_uba_fasta_db(source_file):
    command = ['./scripts/make_uba_fasta_db.sh', source_file.file_name]
    dest_file= source_file.file_name+'.fasta'
    run_msg='making fasta file db for ubas'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def make_fasta_from_blast(blast_result):
    command = ['./scripts/blast_to_fasta.sh', blast_result.file_name]
    dest_file= blast_result.file_name+'.fasta'
    run_msg='making fasta file db for ubas'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)
    
def cluster_at(file_name, percent):
    dest_file=file_name.file_name+'.clstr'+str(percent)
    command = ['cd-hit', '-i', file_name.file_name, '-o', dest_file,
               '-c', str(percent/100), '-n', '5']
    run_msg='making fasta file db for ubas'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, False)

#def cluster_at(file_name, dest_file, percent):
#    command = ['cd-hit', '-i', file_name, '-o', dest_file,
#               '-c', str(percent/100), '-n', '5']
#    run_msg='making fasta file db for ubas'
#    exists=os.path.isfile(dest_file)
#    return run_conditional(command, run_msg, dest_file, exists, False)

def run_mafft(file_name):
    command = ['mafft', file_name.file_name]
    stds = run_command(command)
    dest_file= file_name.file_name+'.aln'
    run_msg='run mafft'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def run_trimal(file_name):
    command = ['trimal', '-automated1' , '-in', file_name.file_name, '-out', file_name.file_name+'.trim']
    dest_file= file_name.file_name+'.trim'
    run_msg='run trimal'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, False)

def run_fasttree(file_name):
    command = ['fasttree', file_name.file_name]
    dest_file= file_name.file_name+'.fast.tree'
    run_msg='run fasttree'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def run_iqtree(file_name):
    command = ['iqtree', '-s', file_name.file_name]
    run_msg='making fasta file db for ubas'
    exists=os.path.isfile(file_name.file_name+'.treefile')
    return run_conditional(command, run_msg, '', exists, False)

def write_file(string, file_name):
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except Error as e:
            print(e)
    with open(file_name, "w") as text_file:
        text_file.write(string)
        text_file.close()
