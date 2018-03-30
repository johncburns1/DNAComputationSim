# This class loads txt files into usable graph data structures
#
#
# author: Jack Burns
# create date: 11/13/2017
# version: 1.0

import networkx


class Graph:

    def __init__(self):
        self.file_path_cities = 'Networks/cities.txt'
        self.file_path_nums = 'Networks/nums.txt'
        self.file_path_peterson = 'Networks/peterson.txt'
        self.citiesGraph = networkx.read_edgelist(self.file_path_cities,
                                                  create_using=networkx.DiGraph(),
                                                  nodetype=str,
                                                  data=[('to', str)])
        self.numsGraph = networkx.read_edgelist(self.file_path_nums,
                                                create_using=networkx.DiGraph(),
                                                nodetype=str,
                                                data=[('to', str)])
        self.petersonGraph = networkx.read_edgelist(self.file_path_peterson,
                                                    create_using=networkx.DiGraph(),
                                                    nodetype=str,
                                                    data=[('to', str)])
