import functions
################################################################################
test_dir = '/home/jax/focal/data/analyses/metagenome/raw/'

#out_group_file = '/home/jax/focal/data/analyses/metagenome/raw/root.fasta'
ref_sequence_path = test_dir + 'canonical/ndm.fasta'
prev_path = test_dir \
            + 'prevalence/protein_fasta_protein_homolog_model_variants.fasta'
output_dir='../output/ndm/'
final_name = 'canon_prev_nrdb_uba.fasta'
def run_all():

    canon_fasta = functions.convert_headers(ref_sequence_path,
                                            output_dir + 'canon/results.fasta',
                                            'canon')
    map_file = functions.make_index_map(ref_sequence_path,
                       output_dir + 'canon/results.map')

    blast_db=functions.make_blastp_db(prev_path,
                       output_dir + 'prevalence/db/prev')
    blast_result = functions.blastp(ref_sequence_path,
            output_dir + 'prevalence/db/prev',
               output_dir + 'prevalence/results.blast',
                                    '1e-80','40')
    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'prev')
    clstr = functions.cluster_at(fasta, 95)
    prev_fasta=clstr
    fasttree = functions.run_fasttree(functions
                                      .run_trimal(functions.run_mafft(clstr)))

    blast_result = functions.blastp_multi(ref_sequence_path,
                     test_dir + 'ncbi-nrdb/',
                     output_dir + 'nrdb/results.blast')
    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'nrdb')
    nrdb_fasta=fasta
    clstr = functions.cluster_at(fasta, 70)
    print(clstr.file_name)
    aln = functions.run_mafft(clstr)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)

    db_fasta_file=functions.make_uba_fasta(test_dir + 'rgi-uba',
                                 'NDM beta-lactamase',
                                 output_dir + 'uba/db.fasta')
    blastdb=functions.make_blastp_db(db_fasta_file.file_name,
                       output_dir + 'uba/db/uba')

    blast_result = functions.blastp(ref_sequence_path,
                               blastdb.file_name,
                                    output_dir + 'uba/result.blast',
                                    '1e-10','40')

    map_file = functions.get_max_info_blast_rep(blast_result)
    fasta = functions.blast_map_to_fasta(map_file, 'uba')
    uba_fasta=fasta
    aln = functions.run_mafft(fasta)
    trim = functions.run_trimal(aln)
    fasttree = functions.run_fasttree(trim)

    #root =  functions.to_result(out_group_file)
    files = [canon_fasta,
                     prev_fasta,
                     nrdb_fasta,
                     uba_fasta]


    final = functions.concatenate_files(files, output_dir + '/' + final_name)

    aln = functions.run_mafft(final)
    trim = functions.run_trimal(aln)
    #fasttree = functions.run_fasttree(trim)
    iqtree = functions.run_iqtree(trim, 'WAG')

    functions.create_table(output_dir
                           + '/sqlite_metatdata', 'ndm_canon')
    functions.insert_all_file(output_dir
                              + '/sqlite_metatdata', 'ndm_canon', output_dir
                              + '/canon/results.map')
    functions.create_table(output_dir
                           + '/sqlite_metatdata', 'ndm_prev')
    functions.insert_all_file(output_dir
                              + '/sqlite_metatdata', 'ndm_prev', output_dir
                              + '/prevalence/results.blast.map')
    functions.make_prev_db(output_dir
                           + '/sqlite_metatdata', 'prev_meta')
    functions.create_table(output_dir
                           + '/sqlite_metatdata', 'ndm_nrdb')
    functions.insert_all_file(output_dir
                              + '/sqlite_metatdata', 'ndm_nrdb', output_dir
                              + '/nrdb/results.blast.map')
    functions.create_table(output_dir
                           + '/sqlite_metatdata', 'ndm_uba')
    functions.insert_all_file(output_dir
                              + '/sqlite_metatdata', 'ndm_uba', output_dir
                              + '/uba/result.blast.map')


    tree=functions.ete_newick_to_tree(functions.car(functions.read_file(iqtree)))

    #tree=functions.read_file('/home/jax/focal/oper/draft/amr_uba/)
    tree2 = functions.annotate(tree, 'ndm' ,output_dir
                               + '/sqlite_metatdata')

    #print(functions.retrieve_gi_info('491516048').std_out)
    functions.render(tree2, output_dir + '/'+final_name+'.svg')

run_all()