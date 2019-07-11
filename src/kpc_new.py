import unittest
import functions

#from ete3 import NCBITaxa

#ncbi = NCBITaxa()
test_dir = '/home/jax/focal/data/analyses/metagenome/raw/'
#canonical  kpc.fasta  kpc.fasta  kpc.fasta  oxa.fasta
# ncbi-nrdb
# prevalence  card_prevalence.txt  protein_fasta_protein_homolog_model_variants.fasta
# rgi-uba
ref_sequence_path = test_dir + 'canonical/kpc.fasta'
prev_path = test_dir + 'prevalence/protein_fasta_protein_homolog_model_variants.fasta'

output_dir='./output_dir/'

def run_all():
    canon_fasta = functions.convert_headers(ref_sequence_path,
                                            output_dir + 'kpc/canon/results.fasta',
                                            'canon')
    map_file = functions.make_index_map(ref_sequence_path,
                       output_dir + 'kpc/canon/results.map')

    blast_db=functions.make_blastp_db(prev_path,
                       output_dir + 'kpc/prevalence/db/prev')
    blast_result = functions.blastp(ref_sequence_path,
            output_dir + 'kpc/prevalence/db/prev',
               output_dir + 'kpc/prevalence/results.blast', '1e-80', '40')
    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'prev')
    clstr = functions.cluster_at(fasta, 99)
    prev_fasta=clstr
    fasttree = functions.run_fasttree(functions.run_trimal(functions.run_mafft(clstr)))

    blast_result = functions.blastp_multi(ref_sequence_path,
                     test_dir + 'ncbi-nrdb/',
                     output_dir + 'kpc/nrdb/results.blast')
    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'nrdb')
    clstr = functions.cluster_at(fasta, 99)
    nrdb_fasta=clstr
    aln = functions.run_mafft(clstr)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
        
        
    db_fasta_file=functions.make_uba_fasta(test_dir + 'rgi-uba',
                                 'KPC beta-lactamase',
                                 output_dir + 'kpc/uba/db.fasta')
    blastdb=functions.make_blastp_db(db_fasta_file.file_name,
                       output_dir + 'kpc/uba/db/uba')

    blast_result = functions.blastp(ref_sequence_path,
                               blastdb.file_name,
                               output_dir + 'kpc/uba/result.blast', '1e-10', '40')

    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'uba')
    uba_fasta=fasta
    aln = functions.run_mafft(fasta)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
        
    filenames = [canon_fasta,
                     prev_fasta,
                     nrdb_fasta,
                     uba_fasta]

    final = functions.concatenate_files(filenames, output_dir + '/kpc/canon_prev_nrdb_uba.fasta')
    
    aln = functions.run_mafft(final)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
    
run_all()
#        aln = run_mafft('./tmp/canon_prev_nrdb_uba.fasta')
#        trim = run_trimal(aln)
#        fasttree = run_fasttree(trim)
