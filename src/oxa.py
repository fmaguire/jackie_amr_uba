import functions
################################################################################
test_dir = '/home/jax/focal/data/analyses/metagenome/raw/'

#out_group_file = '/home/jax/focal/data/analyses/metagenome/raw/root.fasta'
ref_sequence_path = test_dir + 'canonical/oxa.fasta'
prev_path = test_dir \
            + 'prevalence/protein_fasta_protein_homolog_model_variants.fasta'
output_dir='../output/oxa_full/'
final_name = 'canon_prev_nrdb_uba.fasta'
def run_all():

    canon_fasta = functions.convert_headers(ref_sequence_path,
                                            output_dir + 'canon/results.fasta',
                                            'canon')
    map_file = functions.make_index_map(ref_sequence_path,
                       output_dir + 'canon/results.map')


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


    aln = functions.run_mafft(final)
    trim = functions.run_trimal(aln)
    #fasttree = functions.run_fasttree(trim)


    functions.create_table(output_dir
                           + '/sqlite_metatdata', 'oxa_canon')
    functions.insert_all_file(output_dir
                              + '/sqlite_metatdata', 'oxa_canon', output_dir
                              + '/canon/results.map')
    functions.create_table(output_dir
                           + '/sqlite_metatdata', 'oxa_uba')
    functions.insert_all_file(output_dir
                              + '/sqlite_metatdata', 'oxa_uba', output_dir
                              + 'oxa/uba/result.blast.map')

    iqtree = functions.run_iqtree(trim, 'MFP')
    tree=functions.ete_newick_to_tree(functions.car(functions.read_file(iqtree)))

    #tree=functions.read_file('/home/jax/focal/oper/draft/amr_uba/)
    tree2 = functions.annotate(tree, 'oxa' ,output_dir
                               + '/sqlite_metatdata')

    #print(functions.retrieve_gi_info('491516048').std_out)
    functions.render(tree2, output_dir + '/'+final_name+'.svg')

run_all()