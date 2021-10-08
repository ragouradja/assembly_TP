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
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
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
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

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
    count = 1
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
    start_node = get_starting_nodes(graph)
    end_node = get_sink_nodes(graph)
    list_contig = get_contigs(graph, start_node, end_node)
    test_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "test.fna"))
    contig = [("TCAGCGAT", 8), ("TCAGCGAA",8), ("ACAGCGAT", 8), ("ACAGCGAA", 8)]

    save_contigs(list_contig, "fichier.txt")
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # # graphe
    # # Plot the graph
    # draw_graph(graph, "graph_image")
    # #Save the graph in file
    # save_graph(graph, "graph")


if __name__ == '__main__':
    main()
