import unittest
from ete3 import NCBITaxa
from functions import *

ncbi = NCBITaxa()
test_dir = './test/raw/'
ref_sequence_path = test_dir + 'mcr-canonical.fasta' 

class blast_test(unittest.TestCase):

    def test_all(self):
        canon_fasta = convert_headers(ref_sequence_path,
                        './tmp/mcr/canon/results.fasta',
                        'canon')
        map_file = make_index_map(test_dir + 'mcr-canonical.fasta',
                       './tmp/mcr/canon/results.map')

        blast_db=make_blastp_db(test_dir + 'prevalence.fasta',
                       './tmp/mcr/prevalence/db/prev')
        blast_result = blastp(test_dir + 'mcr-canonical.fasta',
               './tmp/mcr/prevalence/db/prev',
               './tmp/mcr/prevalence/results.blast' )
        map_file = get_max_info_blast_rep(blast_result)
        fasta = blast_map_to_fasta(map_file, 'prev')
        clstr = cluster_at(fasta, 95)
        prev_fasta=clstr.file_name
        fasttree = run_fasttree(run_trimal(run_mafft(clstr)))
        
        blast_result = blastp_multi(ref_sequence_path,
                     test_dir + 'ncbi-nrdb/',
                     './tmp/mcr/nrdb/results.blast')
        map_file = get_max_info_blast_rep(blast_result)
        fasta = blast_map_to_fasta(map_file, 'nrdb')
        nrdb_fasta=fasta.file_name
        aln = run_mafft(fasta)
        trim = run_trimal(aln)
        fasttree = run_fasttree(trim)
        
        
        db_fasta_file=make_uba_fasta(test_dir + 'rgi-uba',
                                 'MCR phosphoethanolamine',
                                 './tmp/mcr/uba/db.fasta')
        blastdb=make_blastp_db(db_fasta_file.file_name,
                       './tmp/mcr/uba/db/uba')

        blast_result = blastp(ref_sequence_path,
                               blastdb.file_name,
                               './tmp/mcr/uba/result.blast')

        map_file = get_max_info_blast_rep(blast_result)
        fasta = blast_map_to_fasta(map_file, 'uba')
        uba_fasta=fasta.file_name
        aln = run_mafft(fasta)
        trim = run_trimal(aln)
        fasttree = run_fasttree(trim)
        
        filenames = [canon_fasta.file_name,
                     prev_fasta,
                     nrdb_fasta,
                     uba_fasta]

        with open('./tmp/canon_prev_nrdb_uba.fasta', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

#        aln = run_mafft('./tmp/canon_prev_nrdb_uba.fasta')
#        trim = run_trimal(aln.file_name)
#        fasttree = run_fasttree(trim.file_name)
    '''     
 
                        
    def test_canonical(self):
        convert_headers(ref_sequence_path,
                        './tmp/mcr/canon/results.fasta',
                        'canon')
        make_index_map(ref_sequence_path,
                       './tmp/mcr/canon/results.map')
    def test_prevalence(self):
        make_blastp_db(test_dir + 'prevalence.fasta',
                       './tmp/mcr/prevalence/db/prev')
        blast_result = blastp(ref_sequence_path,
               './tmp/mcr/prevalence/db/prev',
               './tmp/mcr/prevalence/results.blast' )
        map_file = get_max_info_blast_rep(blast_result.file_name)
        fasta = blast_map_to_fasta(map_file.file_name, 'prev')
        aln = run_mafft(fasta.file_name)
        trim = run_trimal(aln.file_name)
        fasttree = run_fasttree(trim.file_name)
        
    def test_nrdb(self):
        blast_result = blastp_multi(ref_sequence_path,
                     test_dir + 'ncbi-nrdb/',
                     './tmp/mcr/nrdb/results.blast')
        map_file = get_max_info_blast_rep(blast_result.file_name)
        fasta = blast_map_to_fasta(map_file.file_name, 'nrdb')
        aln = run_mafft(fasta.file_name)
        trim = run_trimal(aln.file_name)
        fasttree = run_fasttree(trim.file_name)

        
    def test_uba(self):
         db_fasta_file=make_uba_fasta(test_dir + 'rgi-uba',
                                 'MCR phosphoethanolamine',
                                 './tmp/mcr/uba/db.fasta')
         make_blastp_db(db_fasta_file.file_name,
                        './tmp/mcr/uba/db/uba')
         blast_result = blastp(ref_sequence_path,
                               './tmp/mcr/uba/db/uba',
                               './tmp/mcr/uba/results.blast' )
         map_file = get_max_info_blast_rep(blast_result.file_name)
         fasta = blast_map_to_fasta(map_file.file_name, 'uba')
         aln = run_mafft(fasta.file_name)
         trim = run_trimal(aln.file_name)
         fasttree = run_fasttree(trim.file_name)        


         #       blast_map_to_fasta(map_file.file_name, 'uba')
        #    def test_prevalence_blast(self):
#        print(blastp_multi(ref_sequence_path , test_dir + 'ncbi-nrdb/'))
        
        #        make_blastp_db('test/raw/prevalence.fasta', 'tmp/prevalence')
 #       blast_results=blastp('test/raw/mcr-canonical.fasta', 'tmp/prevalence')

#        create_table('tmp/test_db', 'prev')
#        populate_table('tmp/test_db', 'prev', blast_results, extract_aro)
#        print(get_fasta('tmp/test_db', 'prev'))


    #    def test_uba_rgi(self):
#         rgi_card_results = get_rgi_card_results('test/raw/rgi-uba',
#                                                 'MCR phosphoethanolamine')
#         fa_res= make_uba_fasta_db(rgi_card_results.file_name)
#         blast = make_blastp_db(fa_res.file_name, 'tmp/uba')
#         res=blastp('test/raw/mcr-canonical.fasta', blast.file_name)
#         fasta=get_max_info_blast_rep(res.file_name)
#         aln=run_mafft(fasta.file_name)
#         print(aln.file_name)
#         trim=run_trimal(aln.file_name)
#         fasttree=run_fasttree(trim.file_name)
#works#         iqtree=run_iqtree(trim.file_name)

#         cluster_at(fasta.file_name, 100)
#write_file(rgi_card_results, 'tmp/uba.fasta')
#        
#        blast_results = blastp('test/raw/mcr-canonical.fasta', 'tmp/uba')
#        create_table('tmp/test_db', 'uba')
#        populate_table('tmp/test_db', 'uba', blast_results, extract_uba)
#        cluster_at('tmp/uba.fasta', 90)
#        write_file(run_mafft('tmp/uba.fasta.clstr90'), 'tmp/uba.fasta.clstr90.aln')
#        run_trimal('tmp/uba.fasta.clstr90.aln')
#        write_file(run_fasttree('tmp/uba.fasta.clstr90.aln.trim'), 'tmp/uba.fasta.clstr90.aln.trim.tree')
#        run_iqtree('tmp/uba.fasta.clstr90.aln.trim')

        
    def test_prevalence_blast(self):
        make_blastp_db('test/raw/prevalence.fasta', 'tmp/prevalence')
        blast_results=blastp('test/raw/mcr-canonical.fasta', 'tmp/prevalence')
        create_table('tmp/test_db', 'prev')
        populate_table('tmp/test_db', 'prev', blast_results, extract_aro)
        print(get_fasta('tmp/test_db', 'prev'))

    def test_blast_nrdb(self):
        create_table('tmp/test_db', 'nrdb')
        blasty = blastp_multi('test/raw/mcr-canonical.fasta', db_file_list_to_db_names(retrieve_files('test/raw/ncbi-nrdb/', 'nr*[0-9][0-9].*')))
        for result in blasty:
            print(result[0].decode('utf-8').splitlines())
            for blast_result in result[0].decode('utf-8').splitlines():
                info = extract_blast_results(blast_result, extract_gi)
                insert_into_table('tmp/test_db', 'nrdb', info)
        #print_table('tmp/test_db', 'nrdb')
        print(get_fasta('tmp/test_db', 'nrdb'))
        write_file(get_fasta('tmp/test_db', 'nrdb'), 'tmp/tmp.fasta')
    '''
