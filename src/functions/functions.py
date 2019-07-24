import subprocess
import os
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
from ete3 import AttrFace
import pandas as pd

################################################################################
CommandResults = collections.namedtuple('CommandResults', ['file_name',
                                                           'std_out',
                                                           'std_err'])
ncbi = NCBITaxa()
if not os.path.isfile('./taxdump.tar.gz'):
    ncbi.update_taxonomy_database()


def extract_tax_id(ncbiOutput):
    ncbiOutput = ncbiOutput.split('[')
    ncbiOutput = ncbiOutput[1].split(']')[0]
    return (ncbiOutput)


def get_common_ancestor(tax_ids):
    t = ncbi.get_topology(tax_ids, intermediate_nodes=True)
    common = t.get_common_ancestor(tax_ids).name
    common_ancestor = ncbi.get_taxid_translator([common]).values()
    common_ancestor = extract_name(common_ancestor)
    return (common_ancestor)


def extract_name(name):
    name = str(name)
    name = name.split('(')[1]
    name = name.split(')')[0]
    name = name.split('[')[1]
    name = name.split(']')[0]
    name = name.replace("'", '')
    return (name)


def feature_to_face(feature_string, node, col, colour):
    s1 = node.write(features=[feature_string])
    s2 = re.findall(r'' + str(feature_string) + '=[\w\s]+', s1)
    if len(s2) != 0:
        genus = s2[0].split('=')[1]

        hola = TextFace(genus)
        # Set some attributes
        hola.margin_top = 5
        hola.margin_right = 5
        hola.margin_left = 5
        hola.margin_bottom = 5
        hola.opacity = 0.7  # from 0 to 1
        #hola.inner_border.width = 1  # 1 pixel border
        #hola.inner_border.type = 1  # dashed line
        #hola.border.width = 1
        hola.background.color = colour
        if s2[0].split('=')[0] is 'extra':
            hola.bold = True
        node.add_face(hola, column=col, position="aligned")
    else:
        genus = ""
        hola = TextFace(genus)
        # Set some attributes
        hola.margin_top = 5
        hola.margin_right = 5
        hola.margin_left = 5
        hola.margin_bottom = 5
        hola.opacity = 0.7  # from 0 to 1
        # hola.inner_border.width = 1  # 1 pixel border
        # hola.inner_border.type = 1  # dashed line
        # hola.border.width = 1
        hola.background.color = colour

        node.add_face(hola, column=col, position="aligned")


def make_prev_db(db_file, name):
    try:
        conn = sqlite3.connect(db_file)
        prev = pd.read_csv("/home/jocelyn/data/analysis/scripts/index-for-model-sequences.txt", sep='\t')
        prev.to_sql(name, conn, if_exists='replace', index=False)
        curs = conn.cursor()
    except Error as e:
        print(e)
    finally:
        conn.close()


def select_prev_id(db_file, table_name, id_num):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        return (curs.execute(
            "SELECT species_name FROM " + table_name + " WHERE prevalence_sequence_id = " + str(id_num)).fetchall())
    except Error as e:
        print(e)
    finally:
        conn.close()
def annotate(tree, determinant_string, db_file):
    nst1 = NodeStyle()
    nst1["hz_line_color"] = "LightSteelBlue"
    nst1["vt_line_color"] = "LightSteelBlue"
    tree.set_style(nst1)


    colors = ['LightSteelBlue', 'DarkSeaGreen', 'Khaki', 'Moccasin', ]
#colors=[ 'LightSteelBlueGainsboro', 'DarkSeaGreen', 'BlanchedAlmond', 'LightPink',]

    for node in tree.traverse():
        #print(node.name)
        #node.add_feature('monoph_info', str(node.name.split('|')[1]))
        nst1 = NodeStyle()
        nst1["hz_line_color"] = "LightSteelBlue"
        nst1["vt_line_color"] = "LightSteelBlue"
        node.set_style(nst1)
        print(node.name.split('_'))
        if node.is_leaf():
            label_arr = node.name.split('_')
            if (label_arr[1]) == 'nrdb':
                N = AttrFace("name", fsize=12)
                N.background.color = colors[0]
                node.add_face(N, 0, position="aligned")
                N.opacity = 0.2  # from 0 to 1

                #nst1 = NodeStyle()
                #nst1["bgcolor"] = "LightSteelBlue"
                #node.set_style(nst1)
                lcl_id = select_id_from_table(db_file, determinant_string
                                              + '_nrdb', label_arr[2])

                for pair in get_nrdb_lineage(lcl_id):
                    node.add_feature(pair[0].replace(" ", "_"), pair[1])
                    #if pair[0] is 'genus':
                        #node.add_face(AttrFace(node.Genes, text_prefix="/"), column=2, position="branch-top")
                        #node.add_feature('monoph_info', str(node.name.split('|')[1]) + pair[1])
                        #print(str(node.name.split('|')[1]) + pair[1])
                node.add_feature('extra', get_extra_info_nrdb(lcl_id)
                                 .replace("-", " ").replace('.','_'))
                feature_to_face('superkingdom', node, 1, colors[0])
                feature_to_face('phylum', node, 2, colors[0])
                feature_to_face('class', node, 3,colors[0])
                feature_to_face('order', node, 4,colors[0])
                feature_to_face('family', node, 5,colors[0])
                feature_to_face('species', node, 6,colors[0])
                feature_to_face('no_rank', node, 7,colors[0])
                feature_to_face('extra', node, 8, colors[0])

            elif (label_arr[1]) == 'uba':
                N = AttrFace("name", fsize=12)
                N.background.color = colors[1]
                node.add_face(N, 0, position="aligned")
                N.opacity = 0.2  # from 0 to 1
                lcl_id = select_id_from_table(db_file, determinant_string
                                              + '_uba', label_arr[2])
                ubaName = (lcl_id.split('\t')[3]).split('|')[0]
                metadata = pd.read_csv("/home/jocelyn/data/UBA_metadata.tsv",
                                       sep='\t')
                metadata = metadata.drop(metadata.columns[0], axis=1)
                metadata = metadata.drop(metadata.columns[0], axis=1)
                metadata.rename(columns={'UBA Genome ID': 'UBAGenomeID'},
                                inplace=True)
                uba = metadata.loc[metadata.UBAGenomeID == ubaName]
                ubaList = uba.values.tolist()
                ubaList = ubaList[0]
                #print(ubaList[6].rsplit(' ', 1)[0])

                result = ubaList[6].rsplit(' ', 1)[0]
                tax_id = ncbi.get_name_translator([result])[result]
                #print(tax_id)
                lineage = ncbi.get_lineage(tax_id[0])
                names = ncbi.get_taxid_translator(lineage)
                rank = ncbi.get_rank(lineage)

                for pair in list(map(lambda x: (rank[x], names[x]), lineage))[2:]:
                    node.add_feature(pair[0].replace(" ", "_"), pair[1])
                feature_to_face('superkingdom', node, 1, colors[1])
                feature_to_face('phylum', node, 2, colors[1])
                feature_to_face('class', node, 3,colors[1])
                feature_to_face('order', node, 4,colors[1])
                feature_to_face('family', node, 5,colors[1])
                feature_to_face('species', node, 6,colors[1])
                feature_to_face('no_rank', node, 7,colors[1])
                feature_to_face('extra', node, 8, colors[1])
                #nst2 = NodeStyle()
                #nst2["bgcolor"] = "Moccasin"
                #node.set_style(nst2)

            elif (label_arr[1]) == 'canon':
                N = AttrFace("name", fsize=12)
                N.background.color = colors[2]
                node.add_face(N, 0, position="aligned")
                N.opacity = 0.2  # from 0 to 1
                lcl_id = select_id_from_table(db_file, determinant_string
                                              + '_canon', label_arr[2])
                for pair in get_canon_lineage(lcl_id):
                    node.add_feature(pair[0].replace(" ", "_"), pair[1])
                node.add_feature('extra', get_extra_info_canon(lcl_id)
                                 .replace("-", " ").replace('.','_'))
                feature_to_face('superkingdom', node, 1, colors[2])
                feature_to_face('phylum', node, 2, colors[2])
                feature_to_face('class', node, 3,colors[2])
                feature_to_face('order', node, 4,colors[2])
                feature_to_face('family', node, 5,colors[2])
                feature_to_face('species', node, 6,colors[2])
                feature_to_face('no_rank', node, 7,colors[2])
                feature_to_face('extra', node, 8, colors[2])

                #nst2 = NodeStyle()
                #nst2["bgcolor"] = "Orange"
                #node.set_style(nst2)
            elif (label_arr[1]) == 'prev':
                N = AttrFace("name", fsize=12)
                N.background.color = colors[3]
                node.add_face(N, 0, position="aligned")
                N.opacity = 0.2  # from 0 to 1
                lcl_id = select_id_from_table(db_file, determinant_string
                                              + '_prev', label_arr[2])
                line_split = lcl_id.split()[3]

                header_line = line_split.split('|')

                id1 = header_line[0]
                id2 = id1.split(':')[1]

                tax_ids = []
                for result in select_prev_id(db_file, 'prev_meta', id2):
                    tax_id = ncbi.get_name_translator([result[0]])[result[0]]

                    tax_ids.append(str(tax_id[0]))
                print(tax_ids)


                try:

                    t = ncbi.get_topology(tax_ids, intermediate_nodes=True)
                    common = t.get_common_ancestor(tax_ids).name
                    lineage = ncbi.get_lineage(common)  # for now just first

                    names = ncbi.get_taxid_translator(lineage)
                    rank = ncbi.get_rank(lineage)
                    for pair in list(map(lambda x: (rank[x], names[x]), lineage))[2:]:
                        node.add_feature(pair[0].replace(" ", "_"), pair[1])
                    #print(get_extra_info_prev(lcl_id))
                    node.add_feature('extra', get_extra_info_prev(lcl_id)
                                     .replace("-", " ").replace('.', '_'))

                    feature_to_face('superkingdom', node, 1, colors[3])
                    feature_to_face('phylum', node, 2, colors[3])
                    feature_to_face('class', node, 3, colors[3])
                    feature_to_face('order', node, 4, colors[3])
                    feature_to_face('family', node, 5, colors[3])
                    feature_to_face('species', node, 6, colors[3])
                    feature_to_face('no_rank', node, 7, colors[3])
                    feature_to_face('extra', node, 8, colors[3])
                except:
                    lineage = ncbi.get_taxid_translator(tax_ids[0]).name
                    rank = ncbi.get_rank(lineage)
                    for pair in list(map(lambda x: (rank[x], names[x]), lineage))[2:]:
                        node.add_feature(pair[0].replace(" ", "_"), pair[1])
                    #print(get_extra_info_prev(lcl_id))
                    node.add_feature('extra', get_extra_info_prev(lcl_id)
                                     .replace("-", " ").replace('.', '_'))

                    feature_to_face('superkingdom', node, 1, colors[3])
                    feature_to_face('phylum', node, 2, colors[3])
                    feature_to_face('class', node, 3, colors[3])
                    feature_to_face('order', node, 4, colors[3])
                    feature_to_face('family', node, 5, colors[3])
                    feature_to_face('species', node, 6, colors[3])
                    feature_to_face('no_rank', node, 7, colors[3])
                    feature_to_face('extra', node, 8, colors[3])
        else:
            print('root')
        #tree.set_outgroup('lcl|root|1')

        #for node in tree.get_descendants():
        #    if not node.is_leaf() and node.support <= 0.7:
        #        print(node.support)
        #        node.delete()
        #print('checking monophyly')
        #print(node.name)
        #monophs = tree.check_monophyly(values=[], target_attr='monoph_info', ignore_missing=True)
        #print(monophs)

    return tree


def render(tree, name):
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    tree_style.draw_guiding_lines = True
    #tree_style.mode  = "c"
    new_tree_str = tree.write(format=9, features=['genus', 'species',
                                                  'no_rank'])
    tree_style.optimal_scale_level = "full"
    tree_style.scale = 1000
    tree.show_branch_support = True

    write_file(new_tree_str, name + '.ete.tree')
    tree.render(name, tree_style=tree_style)


def read_file(source_file):
    with open(source_file.file_name, "r") as myfile:
        data = myfile.readlines()
    return data


def car(array):
    return array[0]


def ete_newick_to_tree(newick_string):
    return Tree(newick_string, format=0)


def get_canon_lineage(string):
    line_split = string.split('|', 2)
    header_line = line_split[2].split('[')
    organism = header_line[1].split(']')[0]
    tax_id = ncbi.get_name_translator([organism])[organism]
    lineage = ncbi.get_lineage(tax_id[0])  # for now just first
    names = ncbi.get_taxid_translator(lineage)
    rank = ncbi.get_rank(lineage)
    return list(map(lambda x: (rank[x], names[x]), lineage))[2:]


def get_nrdb_lineage(string):
    try:
        tax_id = string.replace('\n','').split('\t').pop().split('-')
        #print(tax_id)
        #print(type(tax_id))
        t = ncbi.get_topology(tax_id, intermediate_nodes=True)
        common = t.get_common_ancestor(tax_id).name
        lineage = ncbi.get_lineage(common)  # for now just first
        names = ncbi.get_taxid_translator(lineage)
        rank = ncbi.get_rank(lineage)
    except:
        tax_id = string.replace('\n','').split('\t').pop().split('-')
        #print(tax_id)
        #print(type(tax_id))
        t = ncbi.get_topology(tax_id, intermediate_nodes=True)
        #common = t.get_common_ancestor(tax_id).name
        lineage = ncbi.get_lineage(tax_id[0])  # for now just first
        names = ncbi.get_taxid_translator(lineage)
        rank = ncbi.get_rank(lineage)

    return (list(map(lambda x: (rank[x], names[x]), lineage))[2:])


def select_id_from_table(db_file, table_name, id_num):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        return (curs.execute("SELECT * FROM " + table_name + " WHERE id = "
                             + id_num).fetchall()[0])[1]
    except Error as e:
        print(e)
    finally:
        conn.close()


def print_table(db_file, table_name):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        for row in curs.execute("SELECT * FROM " + table_name):
            print(str(row))
    except Error as e:
        print(e)
    finally:
        conn.close()


def insert_all_file(db_file, table_name, source_file):
    f = open(source_file, "r")
    for x in f:
        pair = tuple(x.split('\t', 1))
        insert_into_table(db_file, table_name, pair)


def create_table(db_file, table_name):
    try:
        conn = sqlite3.connect(db_file)
        curs = conn.cursor()
        create_table_prev = 'CREATE TABLE ' + table_name \
                            + ' (id integer PRIMARY KEY, string text NOT NULL);'
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
        idst = pair[0].replace(';', '-')
        string = pair[1].replace(';', '-')
        insertion = 'INSERT INTO ' + table_name + '(id, string) VALUES("' \
                    + idst + '","' + string + '")'
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
        print(stds[0].decode('utf-8'))
        if write_truth:
            write_file(stds[0].decode('utf-8'), dest_file)
    else:
        print('exists, not doing action: ' + run_msg)
        stds = (b'', b'')
    return CommandResults(dest_file,
                          stds[0].decode('utf-8'),
                          stds[1].decode('utf-8'))

def to_result(str, str2='', str3=''):
    return CommandResults(str,
                          str2,
                          str3)

def concatenate_files(file_list, dest_file):
    file_list_expanded = list(map(lambda x: x.file_name, file_list))
    command = ['./scripts/general/concatenate_files.sh'] + file_list_expanded
    run_msg = 'concatenating'
    exists = os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)


def convert_headers(source_file, dest_file, label):
    command = ['./scripts/canon/convert_headers.sh',
               source_file,
               label]
    run_msg = 'converting headers'
    exists = os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)


def blastp_multi(query, source_file, dest_file):
    command = ['./scripts/nrdb/blast_nr.sh',
               source_file, query]
    run_msg = 'blasting multiple'
    exists = os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)


def make_uba_fasta(directory, search_string, dest_file):
    command = ['./scripts/uba/make_uba_fasta_db_direct.sh',
               directory, search_string]
    run_msg = 'making uba_fasta'
    exists = os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)


def make_index_map(source_file, dest_file):
    command = ['./scripts/canon/make_index_map.sh',
               source_file]
    run_msg = 'making index map'
    exists = os.path.isfile(dest_file)
    return run_conditional(command, run_msg, dest_file, exists, True)


def blast_map_to_fasta(source_file, label):
    command = ['./scripts/general/blast_to_fasta.sh',
               source_file.file_name, label]
    destination_file = source_file.file_name + '.fasta'
    run_msg = 'running blast map to fasta'
    exists = os.path.isfile(destination_file)
    return run_conditional(command, run_msg, destination_file, exists, True)


def make_blastp_db(source_file, out_database_name):
    command = ['makeblastdb', '-in', source_file, '-input_type', 'fasta',
               '-dbtype', 'prot', '-out', out_database_name]
    exists = (os.path.isfile(out_database_name + ".pin") or
              os.path.isfile(out_database_name + ".phr") or
              os.path.isfile(out_database_name + ".psq"))
    run_msg = 'making blast db'
    return run_conditional(command, run_msg, out_database_name, exists, False)


def blastp(query, source_file, dest_file_name, evalue, coverage):
    command = ['blastp', '-query', query, '-db', source_file,
               '-outfmt', '6 sseqid sseq evalue bitscore qcovsq',
               '-evalue', evalue, '-qcov_hsp_perc', coverage]
    run_msg = 'running blast map to fasta'
    exists = os.path.isfile(dest_file_name)
    return run_conditional(command, run_msg, dest_file_name, exists, True)


def get_max_info_blast_rep(blast_result_file):
    command = ['./scripts/general/max_blast_results_with_index.sh',
               blast_result_file.file_name]
    destination_file = blast_result_file.file_name + '.map'
    run_msg = 'running blast map to fasta'
    exists = os.path.isfile(destination_file)
    return run_conditional(command, run_msg, destination_file, exists, True)


def cluster_at(file_name, percent):
    print(inspect.currentframe().f_code.co_name)
    destination_file = file_name.file_name + '.clstr' + str(percent)
    command = ['cd-hit', '-i', file_name.file_name, '-o', destination_file,
               '-c', str(percent / 100), '-n', '5']
    run_msg = 'making fasta file db for ubas'
    exists = os.path.isfile(destination_file)
    return run_conditional(command, run_msg, destination_file, exists, False)


def retrieve_gi_info(gi_num):
    command = ['./scripts/general/get_gi.sh', str(gi_num)]
    run_msg='retrieving info for gi ' + str(gi_num)
    return run_conditional(command, run_msg, '', False, False)


def run_mafft(file_name):
    command = ['mafft', file_name.file_name]
    destination_file = file_name.file_name + '.aln'
    run_msg = 'run mafft'
    exists = os.path.isfile(destination_file)
    return run_conditional(command, run_msg, destination_file, exists, True)


def run_trimal(file_name):
    command = ['trimal', '-automated1', '-in', file_name.file_name, '-out',
               file_name.file_name + '.trim']
    destination_file = file_name.file_name + '.trim'
    run_msg = 'run trimal'
    exists = os.path.isfile(destination_file)
    return run_conditional(command, run_msg, destination_file, exists, False)


def run_fasttree(file_name):
    command = ['fasttree', file_name.file_name]
    destination_file = file_name.file_name + '.fast.tree'
    run_msg = 'run fasttree'
    exists = os.path.isfile(destination_file)
    return run_conditional(command, run_msg, destination_file, exists, True)


def run_iqtree(file_name):
    command = ['iqtree', '-s', file_name.file_name]
    run_msg = 'running iqtree'
    exists = os.path.isfile(file_name.file_name + '.treefile')
    return run_conditional(command, run_msg, file_name.file_name
                           + '.treefile', exists, False)


def write_file(string, file_name):
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except Error as e:
            print(e)
    with open(file_name, "w") as text_file:
        text_file.write(string)
        text_file.close()

def get_extra_info_canon(canon_id):
    return canon_id.split(' ')[0].split('|')[2]

def get_extra_info_prev(prev_id):
    return prev_id.split('\t')[3].split('|')[1].replace('ARO_Name:', '')

def get_extra_info_nrdb(nrdb_id):
    return nrdb_id.split('\t')[2]