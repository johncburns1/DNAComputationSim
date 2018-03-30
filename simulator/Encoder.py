# This is an implementation of the DNA encoding algorithm takes an input network and encodes
# it into a DNA sequence using an Adleman-like algorithm
#
# author: Jack Burns
# create date: 11/13/2017
# version 1.0

from Graph import Graph
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
import random as rand


class Encoder:
    A_KEY = "A"
    T_KEY = "T"
    C_KEY = "C"
    G_KEY = "G"

    SEQ_LEN = 20
    COMP_LEN = SEQ_LEN

    base_list = [A_KEY, T_KEY, C_KEY, G_KEY]

    def encodeNodes(self, graph):
        node_dna_map = {}
        for node in graph.nodes():
            sequence = ""
            for i in range(self.SEQ_LEN):
                sequence += str((self.base_list[rand.randint(0, len(self.base_list) - 1)]))
            node_dna_map[node] = sequence
        return node_dna_map

    def encodeEdges(self, graph, encoded_nodes):
        encoded_edges = {}
        edges = graph.edges()

        for xnode in graph.nodes():
            for ynode in graph.nodes():
                if (xnode, ynode) in edges:
                    s1 = encoded_nodes[xnode]
                    s2 = encoded_nodes[ynode]
                    encoded_edges[(xnode, ynode)] = s1[-10:] + s2[:10]

        return encoded_edges

    def generateComplements(self, encoded_nodes):
        complements = {}
        for node in encoded_nodes:
            comp_sequence = self.getSeqComplement(encoded_nodes[node])
            complements[node] = comp_sequence

        return complements

    def getComplement(self, nuc):
        if nuc == self.A_KEY:
            return self.T_KEY
        elif nuc == self.T_KEY:
            return self.A_KEY
        elif nuc == self.G_KEY:
            return self.C_KEY
        elif nuc == self.C_KEY:
            return self.G_KEY
        else:
            return None

    def getSeqComplement(self, seq):
        str = ""
        for s in seq:
            str = str + self.getComplement(s)
        return str[::-1]

    def toDSEQ(self, graph, edges, nodes):
        complements = self.generateComplements(nodes)
        dna = []
        offset = -10
        for edge in edges:
            seq = Dseq(edges[edge], complements[edge[1]], ovhg=offset)
            x = Dseqrecord(seq)
            x.name = edge[0] + "_" + edge[1]
            x.seq = seq
            dna.append(x)

        return dna
