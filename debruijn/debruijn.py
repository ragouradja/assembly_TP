#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
from decimal import getcontext
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter, le, length_hint, sub
import random
random.seed(9001)
from random import randint
import statistics
import pickle

__author__ = "Ragousandirane Radjasandirane"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Ragousandirane Radjasandirane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ragousandirane Radjasandirane"
__email__ = "radja.ragou@gmail.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as filin:
        content = filin.readlines()
        for i in range(1,len(content),4):
            yield content[i].strip()


def cut_kmer(read, kmer_size):
    for i in range(len(read) - (kmer_size - 1)):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    gen_seq = read_fastq(fastq_file)
    seq = gen_seq
    dict_kmer = {}
    for read in seq:
        kmer_list = cut_kmer(read, kmer_size)
        for kmer in kmer_list:
            if kmer not in dict_kmer:
                dict_kmer[kmer] = 1
            else:
                dict_kmer[kmer] += 1
    return dict_kmer


def build_graph(kmer_dict):
    G = nx.DiGraph()
    for kmer,weight in kmer_dict.items():
        value1 = kmer[:-1]
        value2 = kmer[1:]
        G.add_edge(value1, value2, weight = weight)
    return G

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node:
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    std_weight = std(weight_avg_list)
    if std_weight > 0:
        max_weight = max(weight_avg_list)
        max_index = weight_avg_list.index(max_weight)
        path_list.remove(path_list[max_index])
        remove_paths(graph, path_list,
        delete_entry_node, delete_sink_node)

    elif std_weight == 0:
        std_length = std(path_length)
        if std_length > 0:
            max_length = max(path_length)
            max_index_length = path_length.index(max_length)
            path_list.remove(path_list[max_index_length])
            remove_paths(graph, path_list,
            delete_entry_node, delete_sink_node)

        elif std_length == 0:
            r_index = random.randint(0, len(path_length))
            path_list.remove(path_list[r_index])
            remove_paths(graph, path_list,
            delete_entry_node, delete_sink_node)
    return graph
    
def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    if nx.has_path(graph,ancestor_node, descendant_node):
        length_list = []
        weight_avg_list = []
        path = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
        for subpath in path:
            length_list.append(len(subpath))
            weight_avg_list.append(path_average_weight(graph, subpath))
        if len(weight_avg_list) > 1:
            return select_best_path(graph, path, length_list, weight_avg_list)
    return graph

def simplify_bubbles(graph):
    bubble = False

    for node in graph.nodes:
        node_predecessor = list(graph.predecessors(node))
        nb_pred = len(node_predecessor)
        if nb_pred > 1:
            for i in range(nb_pred - 1):
                for j in range(i + 1, nb_pred):
                    node_ancestor = nx.lowest_common_ancestor(graph, node_predecessor[i], node_predecessor[j])
                    if node_ancestor:
                        bubble = True
                        break
            if bubble:
                graph = simplify_bubbles(solve_bubble(graph, node_ancestor, node))
                break
    return graph

def solve_entry_tips(graph, starting_nodes):
    for node in graph.nodes:
        all_predecessors = list(nx.edge_dfs(graph,node, orientation='reverse'))
        if len(all_predecessors) > 1:
            first_path = [node]
            second_path = [node] 
            first_connected = False
            second_connected = False
            for predecessor in all_predecessors:
                if not first_connected:
                    first_path.append(predecessor[0])
                else:
                    second_path.append(predecessor[0])
                if predecessor[0] in starting_nodes:
                    if not first_connected:
                        first_connected = True
                    else:
                        second_connected = True
            if first_connected and second_connected:
                first_path.reverse()
                second_path.reverse()
                length_list = [len(first_path), len(second_path)]
                weight_list = [path_average_weight(graph, first_path),
                path_average_weight(graph, second_path)]
                graph = select_best_path(graph, [first_path, second_path],
                length_list, weight_list, delete_entry_node = True)
                solve_entry_tips(graph, starting_nodes)
                break
    return graph

def solve_out_tips(graph, ending_nodes):
    for node in list(graph.nodes)[::-1]:
        all_successors = list(nx.edge_dfs(graph,node))
        if len(all_successors) > 1:
            first_path = [node]
            second_path = [node] 
            first_connected = False
            second_connected = False
            for successor in all_successors:
                if not first_connected:
                    first_path.append(successor[1])
                else:
                    second_path.append(successor[1])
                if successor[1] in ending_nodes:
                    if not first_connected:
                        first_connected = True
                    else:
                        second_connected = True
            if first_connected and second_connected:
                length_list = [len(first_path), len(second_path)]
                weight_list = [path_average_weight(graph, first_path),
                path_average_weight(graph, second_path)]
                graph = select_best_path(graph, [first_path, second_path],
                length_list, weight_list, delete_sink_node = True)
                solve_entry_tips(graph, ending_nodes)
                break
    return graph


def get_starting_nodes(graph):
    entry_node = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            entry_node.append(node)
    return entry_node

def get_sink_nodes(graph):
    sink_node = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            sink_node.append(node)
    return sink_node

def get_contigs(graph, starting_nodes, ending_nodes):
    contig_list = []
    for start in starting_nodes:
        for end in ending_nodes:
           if nx.has_path(graph,start, end):
               path = list(nx.all_simple_paths(graph, start, end))[0]
               base = path[0]
               seq = base
               len_path = len(path)
               for i in range(1,len_path):
                   seq += path[i][-1]
               final_len = len(seq)
               contig_list.append((seq, final_len))
    return contig_list

def save_contigs(contigs_list, output_file):
    count = 0
    with open(output_file, "wt") as filout:  
        for contigs in contigs_list:
            path = contigs[0]
            len_contig = contigs[1]
            fasta_seq = fill(path)
            header = f">contig_{count} len={len_contig}"
            filout.write(f"{header}\n")
            filout.write(f"{fasta_seq}\n")
            count += 1

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    dico = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dico)
 
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    graph = simplify_bubbles(graph)
    graph_without_entry_tips = solve_entry_tips(graph, starting_nodes)
    graph_final = solve_out_tips(graph_without_entry_tips, sink_nodes)


    list_contigs = get_contigs(graph_final, starting_nodes, sink_nodes)
    save_contigs(list_contigs, args.output_file)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # # graphe
    # # Plot the graph
    # draw_graph(graph, "graph_image")
    # #Save the graph in file
    # save_graph(graph, "graph")


if __name__ == '__main__':
    main()
