# This is the main method to run the DNA encoding algorithm
#
#
# author: Jack Burns
# create date: 11/13/2017
# version 1.0
import sys

import networkx
from pydna.assembly import Assembly
from Encoder import Encoder
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import Anneal
import matplotlib.pyplot as plt
import networkx as nx

enc = Encoder()


def main():
    if len(sys.argv) < 4:
        print("Input is not valid")
        return exit(1)
    else:
        file_prefix = 'Networks/'
        if file_prefix not in str(sys.argv[1]):
            file_path = 'Networks/' + str(sys.argv[1])
        else:
            file_path = str(sys.argv[1])
        try:
            graph = networkx.read_edgelist(file_path,
                                           create_using=networkx.DiGraph(),
                                           nodetype=str,
                                           data=[('to', str)])
        except FileNotFoundError:
            print("Wrong file or file path")
            return exit(1)

        else:
            start_vertex = str(sys.argv[2]).upper()
            end_vertex = str(sys.argv[3]).upper()
            nodes = enc.encodeNodes(graph)

            if start_vertex not in nodes:
                print("could not identify start vertex")
                return exit(1)
            elif end_vertex not in nodes:
                print("could not identify end vertex")
                return exit(1)
            else:
                ham_graph_path = 'graphs/ham_path'
                edges = enc.encodeEdges(graph, nodes)
                node_names = getFilterNodeList(nodes, start_vertex, end_vertex)
                nx.draw_circular(graph, with_labels=True)
                plt.savefig('graphs/graph.png')
                plt.clf()

                x = AssembleNAnneal(graph, nodes, edges, start_vertex, end_vertex)
                filtered_paths = filter(x, nodes, node_names)
                results = getResults(nodes, edges, filtered_paths, ham_graph_path)

                if results[0] == True:
                    print("\nTHERE IS A HAMILTONIAN PATH " + str(results[1]))
                    return exit(0)

                else:
                    print("\nTHERE IS NOT HAMILTONIAN PATH")
                    return exit(0)


def getResults(nodes, edges, filtered_results, file_name):
    ham_path = False
    ham_path_tup = (0, 0)
    if len(filtered_results) != 0:
        for result in filtered_results:
            print("\nHAMILTONIAN PATH: " + str(result.seq))
            ham_path_tup = extractEdges(str(result.seq), nodes, edges)
            nx.draw_circular(ham_path_tup[0], with_labels=True)
            plt.savefig(file_name)
            plt.clf()
        ham_path = True
    return (ham_path, ham_path_tup[1])


def AssembleNAnneal(graph, nodes, edges, start, end):
    dseq_list = enc.toDSEQ(graph, edges, nodes)
    p1 = Dseqrecord(nodes[start])
    p2 = Dseqrecord(enc.getSeqComplement(nodes[end]))
    assembly = Assembly(dseq_list, limit=10, only_terminal_overlaps=True)
    print("\n" + str(assembly) + "\n")
    candidates = []

    for i in range(len(assembly.linear_products)):
        product = assembly.linear_products[i]
        template = Dseqrecord(product)
        pcr = Anneal([p1, p2], template, limit=10)
        gel = len(nodes) * enc.SEQ_LEN

        if len(pcr.products) != 0:
            print(product.detailed_figure())
            print(product.figure())
            for p in pcr.products:
                if len(p.seq) == gel:
                    p.seq = p.seq[10:]
                    p.seq = p.seq[:-10]
                    candidates.append(p)

    # print("\n" +str(nodes))
    # print(str(edges) +"\n")
    return candidates


def filter(paths, nodes, node_name_list):
    n = len(node_name_list)
    unique_paths = list(set(paths))
    gp = []
    pp = {}

    for path in unique_paths:
        pp[path] = 0

    while len(node_name_list) > 0:
        x = node_name_list.pop()
        for path in unique_paths:
            if nodes[x] in path:
                pp[path] += 1

    for path in pp:
        if pp[path] == n:
            gp.append(path)

    return gp


def extractEdges(path, node_list, edge_list):
    g = nx.DiGraph()
    order_edges = [0] * int(len(path) / enc.SEQ_LEN)
    for node in node_list:
        g.add_node(node)

    for edge in edge_list:
        if edge_list[edge] in path:
            index = int(path.index(edge_list[edge]) / enc.SEQ_LEN)
            order_edges[index] = edge
            g.add_edge(edge[0], edge[1])

    return (g, order_edges)


def getFilterNodeList(nodes, start, end):
    node_list = []
    for node in nodes:
        if node != start and node != end:
            node_list.append(node)

    return node_list


if __name__ == '__main__':
    main()
