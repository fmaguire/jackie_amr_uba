import unittest
import functions

#from ete3 import NCBITaxa

#ncbi = NCBITaxa()
test_dir = '/home/jax/focal/data/analyses/metagenome/raw/'
#canonical  kpc.fasta  mcr.fasta  ndm.fasta  oxa.fasta
# ncbi-nrdb
# prevalence  card_prevalence.txt  protein_fasta_protein_homolog_model_variants.fasta
# rgi-uba
ref_sequence_path = test_dir + 'canonical/oxa.fasta'
prev_path = test_dir + 'prevalence/protein_fasta_protein_homolog_model_variants.fasta'

output_dir='./output_dir/'

def run_all():
    canon_fasta = functions.convert_headers(ref_sequence_path,
                                            output_dir + 'oxa/canon/results.fasta',
                                            'canon')
    clstr_canon = functions.cluster_at(canon_fasta, 99)
    map_file = functions.make_index_map(ref_sequence_path,
                       output_dir + 'oxa/canon/results.map')


        
    db_fasta_file=functions.make_uba_fasta(test_dir + 'rgi-uba',
                                 'OXA beta-lactamase',
                                 output_dir + 'oxa/uba/db.fasta')
    blastdb=functions.make_blastp_db(db_fasta_file.file_name,
                       output_dir + 'oxa/uba/db/uba')

    blast_result = functions.blastp(ref_sequence_path,
                               blastdb.file_name,
                               output_dir + 'oxa/uba/result.blast', '1e-180', '90')

    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'uba')
    clstr = functions.cluster_at(fasta, 97)
    uba_fasta=clstr
    aln = functions.run_mafft(clstr)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
        
    filenames = [clstr_canon,
                     uba_fasta]

    final = functions.concatenate_files(filenames, output_dir + '/oxa/canon_uba.fasta')

    aln = functions.run_mafft(final)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
run_all()
#        aln = run_mafft('./tmp/canon_prev_nrdb_uba.fasta')
#        trim = run_trimal(aln.file_name)
#        fasttree = run_fasttree(trim.file_name)
