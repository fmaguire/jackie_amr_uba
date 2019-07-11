import subprocess
import os
import signal
import glob
import re
import sqlite3                                                                                                                                         
from sqlite3 import Error
import collections 
import inspect
from ete3 import NCBITaxa
from ete3 import Tree
from ete3 import TreeStyle
from ete3 import NodeStyle
from ete3 import TextFace
import pandas as pd

CommandResults = collections.namedtuple('CommandResults',['file_name',
                                                          'std_out',
                                                          'std_err'])

ncbi = NCBITaxa()

def get_prev_stuff(aroString,db):
    prev = pd.read_csv("/home/jocelyn/data/analysis/scripts/index-for-model-sequences.txt", sep='\t')
    prev.to_sql('prev_stuff', db, if_exists='append', index=False)
    prevSubset = prev.loc[prev["aro_accession"] == aroString]
    prevSubset = prevSubset.loc[prevSubset["rgi_criteria"] == "Perfect"]
    species = prevSubset["species_name"].unique()
    #print(species)
    taxids = []
    for i in species:
        t = str(getTaxid(i))
        lineage = ncbi.get_lineage(t)
        #print(ncbi.get_rank(lineage))
        taxids.append(t)
    ancestor = getCommonAncestor(taxids)
    return ancestor

def getTaxid(organism):
    taxid = str(ncbi.get_name_translator([organism])[organism])
    taxid = extractTaxid(taxid)
    return(taxid)

def extractTaxid(ncbiOutput):
    ncbiOutput = ncbiOutput.split('[')
    ncbiOutput = ncbiOutput[1].split(']')[0]
    return(ncbiOutput)

def getCommonAncestor(taxids):
    t = ncbi.get_topology(taxids, intermediate_nodes=True)
    common = t.get_common_ancestor(taxids).name
    commonAncestor = ncbi.get_taxid_translator([common]).values()
    commonAncestor = extractName(commonAncestor)
    return(commonAncestor)

def extractName(name):
    name = str(name)
    name = name.split('(')[1]
    name = name.split(')')[0]
    name = name.split('[')[1]
    name = name.split(']')[0]
    name = name.replace("'", '')
    return(name)

def getLineage(organism):
    taxid = str(ncbi.get_name_translator([organism])[organism])
    taxid = extractTaxid(taxid)
    lineage = ncbi.get_lineage(taxid)
    nameOfLineage = ncbi.get_taxid_translator(lineage)
    namesOfLineage = extractLineageNames(nameOfLineage.values())
    return(namesOfLineage)

def extractLineageNames(lineage):
    lineage = str(lineage)
    lineage = lineage.split('(')[1]
    lineage = lineage.split(')')[0]
    lineage = lineage.split('[')[1]
    lineage = lineage.split(']')[0]
    lineage = lineage.replace("'", '')
    return(lineage)

def feature_to_face(feature_string, node, col):
    s1=node.write(features=[feature_string])
    s2=re.findall(r''+str(feature_string)+'=[\w\s]+', s1)
    #print(len(s2))
    if len(s2) != 0:
        genus=s2[0].split('=')[1]

        hola = TextFace(genus)
        # Set some attributes
        hola.margin_top = 10
        hola.margin_right = 10
        hola.margin_left = 10
        hola.margin_bottom = 10
        hola.opacity = 0.5 # from 0 to 1
        hola.inner_border.width = 1 # 1 pixel border
        hola.inner_border.type = 1  # dashed line
        hola.border.width = 1
        hola.background.color = "LightGreen"

        node.add_face(hola, column=col, position = "aligned")


def make_prev_db(db_file, name):

    try:
        conn = sqlite3.connect(db_file)
        prev = pd.read_csv("/home/jocelyn/data/analysis/scripts/index-for-model-sequences.txt", sep='\t')
        prev.to_sql(name, conn, if_exists='replace', index=False)
        curs = conn.cursor()
        #return (curs.execute("SELECT * FROM "+table_name + " WHERE id = " + id_num).fetchall()[0])[1]
    except Error as e:
        print(e)
    finally:
        conn.close()

def select_prev_id(db_file, table_name, id_num):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        return (curs.execute("SELECT species_name FROM "+table_name + " WHERE prevalence_sequence_id = " + str(id_num)).fetchall())
    except Error as e:
        print(e)
    finally:
        conn.close()

def annotate(tree, determinant_string , db_file):

    for node in tree.traverse():
        if node.is_leaf():
            label_arr = node.name.split('|')
           # print(label_arr)
            if(label_arr[1]) == 'nrdb':
                nst1 = NodeStyle()
                nst1["bgcolor"] = "LightSteelBlue"
                node.set_style(nst1)
                #fix this later
                lcl_id = select_id_from_table(db_file, determinant_string + '_nrdb', label_arr[2])

                for pair in get_nrdb_lineage(lcl_id):
                    node.add_feature(pair[0].replace(" ", "_"),pair[1])
                feature_to_face('superkingdom', node, 1)
                feature_to_face('phylum', node, 2)
                feature_to_face('class', node, 3)
                feature_to_face('order', node, 4)
                feature_to_face('family', node, 5)
                feature_to_face('species', node, 6)
                feature_to_face('no_rank', node, 7)


            if(label_arr[1]) == 'uba':
                lcl_id = select_id_from_table(db_file, determinant_string + '_uba', label_arr[2])
                ubaName=(lcl_id.split('\t')[3]).split('|')[0]
                #print(lcl_id)
                #ubaName = "UBA4"
                metadata = pd.read_csv("/home/jocelyn/data/UBA_metadata.tsv", sep='\t')
                metadata = metadata.drop(metadata.columns[0], axis=1)
                metadata = metadata.drop(metadata.columns[0], axis=1)
                metadata.rename(columns = {'UBA Genome ID': 'UBAGenomeID'}, inplace=True)
                uba = metadata.loc[metadata.UBAGenomeID == ubaName]
                ubaList = uba.values.tolist()
                ubaList = ubaList[0]
                print(ubaList[6].rsplit(' ', 1)[0])

                result = ubaList[6].rsplit(' ', 1)[0]
                taxid = ncbi.get_name_translator([result])[result]
                print(taxid)
                lineage = ncbi.get_lineage(taxid[0])
                names = ncbi.get_taxid_translator(lineage)
                rank = ncbi.get_rank(lineage)

                for pair in list(map(lambda x: (rank[x], names[x]), lineage))[2:]:
                    node.add_feature(pair[0].replace(" ", "_"), pair[1])
                feature_to_face('superkingdom', node, 1)
                feature_to_face('phylum', node, 2)
                feature_to_face('class', node, 3)
                feature_to_face('order', node, 4)
                feature_to_face('family', node, 5)
                feature_to_face('species', node, 6)
                feature_to_face('no_rank', node, 7)
                nst2 = NodeStyle()
                nst2["bgcolor"] = "Moccasin"
                node.set_style(nst2)

            if (label_arr[1]) == 'canon':
                lcl_id = select_id_from_table(db_file, determinant_string + '_canon', label_arr[2])
                for pair in get_canon_lineage(lcl_id):
                    node.add_feature(pair[0].replace(" ", "_"),pair[1])
                feature_to_face('superkingdom', node, 1)
                feature_to_face('phylum', node, 2)
                feature_to_face('class', node, 3)
                feature_to_face('order', node, 4)
                feature_to_face('family', node, 5)
                feature_to_face('species', node, 6)
                feature_to_face('no_rank', node, 7)

                nst2 = NodeStyle()
                nst2["bgcolor"] = "Orange"
                node.set_style(nst2)
            if (label_arr[1]) == 'prev':
                lcl_id = select_id_from_table(db_file, determinant_string + '_prev', label_arr[2])
                #id=get_prev_lineage(lcl_id)
                # print(string)
                lineSplit = lcl_id.split()[3]
                # print(lineSplit)
                headerLine = lineSplit.split('|')
                # canonicalSeq = headerLine[0]
                id1 = headerLine[0]
                id2= id1.split(':')[1]
                # print((list(map(lambda x: (rank[x], names[x]), lineage))[2:]))
                #print(id1.split(':')[1])
                taxids=[]
                for result in select_prev_id(db_file, 'prev_meta', id2):
                    #print(result[0])
                    taxid = ncbi.get_name_translator([result[0]])[result[0]]
                    #print(taxid)
                    taxids.append(str(taxid[0]))
                #print(taxids)
                t = ncbi.get_topology(taxids, intermediate_nodes=True)
                common = t.get_common_ancestor(taxids).name

                lineage = ncbi.get_lineage(common) # for now just first

                names = ncbi.get_taxid_translator(lineage)
                rank = ncbi.get_rank(lineage)

                for pair in list(map(lambda x: (rank[x], names[x]), lineage))[2:]:
                    node.add_feature(pair[0].replace(" ", "_"),pair[1])
                feature_to_face('superkingdom', node, 1)
                feature_to_face('phylum', node, 2)
                feature_to_face('class', node, 3)
                feature_to_face('order', node, 4)
                feature_to_face('family', node, 5)
                feature_to_face('species', node, 6)
                feature_to_face('no_rank', node, 7)

    return tree

def render(tree):
    tree_style = TreeStyle()

    tree_style.show_leaf_name = True
    new_tree_str = tree.write(format=9, features=['genus','species', 'no_rank'])
   # print(new_tree_str)
    #new_tree=Tree(new_tree_str)
    write_file(new_tree_str, './mytree.tree')

    tree.render("mytree2.svg", tree_style=tree_style)

def read_file(source_file):
    with open (source_file.file_name, "r") as myfile:
        data=myfile.readlines()
    return data

def car(array):
    return array[0]

def ete_treeify(newick_string):
    return Tree(newick_string, format=0)

if not os.path.isfile('./taxdump.tar.gz'):
    ncbi.update_taxonomy_database()

def get_canon_lineage(string):

    lineSplit = string.split('|', 2)
    headerLine = lineSplit[2].split('[')
    #canonicalSeq = headerLine[0]
    organism = headerLine[1].split(']')[0]
    #print(organism)
    taxid = ncbi.get_name_translator([organism])[organism]

    lineage = ncbi.get_lineage(taxid[0]) # for now just first

    names = ncbi.get_taxid_translator(lineage)
    rank = ncbi.get_rank(lineage)
    #print((list(map(lambda x: (rank[x], names[x]), lineage))[2:]))

    return (list(map(lambda x: (rank[x], names[x]), lineage))[2:])#list(names.values())[1:]


def get_prev_lineage(string):
    #print(string)
    lineSplit = string.split()[3]
    #print(lineSplit)
    headerLine = lineSplit.split('|')
    #canonicalSeq = headerLine[0]
    id = headerLine[0]
    #print((list(map(lambda x: (rank[x], names[x]), lineage))[2:]))
    print(id.split(':')[1])
    #print(select_prev_id(db_file, 'prev_meta', id))
    #print(get_prev_stuff(aro))
    #return (list(map(lambda x: (rank[x], names[x]), lineage))[2:])#list(names.values())[1:]

def get_uba_lineage(string):
    uba=string.split
#    print(uba)
    return []

def get_nrdb_lineage(string):
    taxid=string.split('\t').pop().split('-')
  #  print(taxid[0])
    lineage = ncbi.get_lineage(taxid[0]) # for now just first


    names = ncbi.get_taxid_translator(lineage)
    rank = ncbi.get_rank(lineage)
    return (list(map(lambda x: (rank[x], names[x]), lineage))[2:])#list(names.values())[1:]

def select_id_from_table(db_file, table_name, id_num):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        return (curs.execute("SELECT * FROM "+table_name + " WHERE id = " + id_num).fetchall()[0])[1]
    except Error as e:
        print(e)
    finally:
        conn.close()

def print_table(db_file, table_name):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        for row in curs.execute("SELECT * FROM "+table_name):
            print(str(row))
    except Error as e:
        print(e)
    finally:
        conn.close()

def insert_all_file(db_file, table_name, source_file):
    f = open(source_file, "r")
    for x in f:
        pair = tuple(x.split('\t',1))
        insert_into_table(db_file, table_name, pair)

def create_table(db_file, table_name):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        create_table_prev = 'CREATE TABLE ' + table_name + ' (id integer PRIMARY KEY, string text NOT NULL);'
        #curs.execute('DROP TABLE IF EXISTS ' + table_name )
        curs.execute(create_table_prev)
        conn.commit()
    except Error as e:
        print(e)
    finally:
        conn.close()

def insert_into_table(db_file, table_name, pair):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        idst = pair[0].replace(';','-')
        string = pair[1].replace(';','-')
        insertion = 'INSERT INTO ' + table_name + '(id, string) VALUES("' + idst + '","' + string + '")'
        curs.execute(insertion)
        conn.commit()
    except Error as e:
        print(e)
    finally:
        conn.close()
        
def run_command(command):
    return subprocess.Popen(command,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

def run_conditional(command, run_msg, dest_file, run_truth, write_truth):
    if not run_truth:
        print(run_msg)
        stds = run_command(command)
        #print(stds[0].decode('utf-8'))
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
    print(inspect.currentframe().f_code.co_name)
    command = ['./scripts/general/make_uba_fasta_db.sh', source_file.file_name]
    dest_file= source_file.file_name+'.fasta'
    run_msg='making fasta file db for ubas'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)

def make_fasta_from_blast(blast_result):
    print(inspect.currentframe().f_code.co_name)
    command = ['./scripts/blast_to_fasta.sh', blast_result.file_name]
    dest_file= blast_result.file_name+'.fasta'
    run_msg='making fasta file db for ubas'
    exists=os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)
    
def cluster_at(file_name, percent):
    print(inspect.currentframe().f_code.co_name)
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
    run_msg='running iqtree'
    exists=os.path.isfile(file_name.file_name+'.treefile')
    return run_conditional(command, run_msg, file_name.file_name + '.treefile', exists, False)

def write_file(string, file_name):
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except Error as e:
            print(e)
    with open(file_name, "w") as text_file:
        text_file.write(string)
        text_file.close()
