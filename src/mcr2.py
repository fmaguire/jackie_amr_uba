import functions

test_dir = '/home/jax/focal/data/analyses/metagenome/raw/'

ref_sequence_path = test_dir + 'canonical/mcr.fasta'
prev_path = test_dir + 'prevalence/protein_fasta_protein_homolog_model_variants.fasta'

output_dir='../mcr2_output/'

def run_all():
    canon_fasta = functions.convert_headers(ref_sequence_path,
                                            output_dir + 'mcr2/canon/results.fasta',
                                            'canon')
    map_file = functions.make_index_map(ref_sequence_path,
                       output_dir + 'mcr2/canon/results.map')

    blast_db=functions.make_blastp_db(prev_path,
                       output_dir + 'mcr2/prevalence/db/prev')
    blast_result = functions.blastp(ref_sequence_path,
            output_dir + 'mcr2/prevalence/db/prev',
               output_dir + 'mcr2/prevalence/results.blast',
                                    '1e-180','90')
    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'prev')
    clstr = functions.cluster_at(fasta, 95)
    prev_fasta=clstr
    fasttree = functions.run_fasttree(functions.run_trimal(functions.run_mafft(clstr)))

    blast_result = functions.blastp_multi(ref_sequence_path,
                     test_dir + 'ncbi-nrdb/',
                     output_dir + 'mcr2/nrdb/results.blast')
    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'nrdb')
    nrdb_fasta=fasta
    clstr = functions.cluster_at(fasta, 70)
    print(clstr.file_name)
    aln = functions.run_mafft(clstr)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)

    db_fasta_file=functions.make_uba_fasta(test_dir + 'rgi-uba',
                                 'MCR phosphoethanolamine',
                                 output_dir + 'mcr/uba/db.fasta')
    blastdb=functions.make_blastp_db(db_fasta_file.file_name,
                       output_dir + 'mcr2/uba/db/uba')

    blast_result = functions.blastp(ref_sequence_path,
                               blastdb.file_name,
                                    output_dir + 'mcr2/uba/result.blast',
                                    '1e-180','90')

    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'uba')
    uba_fasta=fasta
    aln = functions.run_mafft(fasta)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
        
    files = [canon_fasta,
                     prev_fasta,
                     nrdb_fasta,
                     uba_fasta]

    final = functions.concatenate_files(files, output_dir + '/mcr2/canon_prev_nrdb_uba.fasta')

    aln = functions.run_mafft(final)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)
    iqtree = functions.run_iqtree(trim)

    functions.create_table(output_dir + '/sqlite_metatdata', 'mcr_canon')
    functions.insert_all_file(output_dir + '/sqlite_metatdata', 'mcr_canon', output_dir + '/mcr2/canon/results.map')
    functions.create_table(output_dir + '/sqlite_metatdata', 'mcr_prev')
    functions.insert_all_file(output_dir + '/sqlite_metatdata', 'mcr_prev', output_dir + '/mcr2/prevalence/results.blast.map')
    functions.make_prev_db(output_dir + '/sqlite_metatdata', 'prev_meta')
    functions.create_table(output_dir + '/sqlite_metatdata', 'mcr_nrdb')
    functions.insert_all_file(output_dir + '/sqlite_metatdata', 'mcr_nrdb', output_dir + '/mcr2/nrdb/results.blast.map')
    functions.create_table(output_dir + '/sqlite_metatdata', 'mcr_uba')
    functions.insert_all_file(output_dir + '/sqlite_metatdata', 'mcr_uba', output_dir + '/mcr2/uba/result.blast.map')

    tree=functions.ete_newick_to_tree(functions.car(functions.read_file(fasttree)))

    tree2 = functions.annotate(tree, 'mcr' ,output_dir + '/sqlite_metatdata')
    
    functions.render(tree2, output_dir + 'mcr2.svg')

run_all()
