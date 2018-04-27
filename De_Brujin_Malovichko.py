from Bio import SeqIO
import Bio
from graphviz import Digraph
import argparse
import pydot

class Vertex:

    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}

    def increase_coverage(self):
        self.coverage += 1


class Edge:

    def __init__(self, k1, k2):
        self.seq = k1 + k2[-1]
        self.n = 2
        self.coverage = 0

    def calc_coverage(self, c1, c2):
        self.coverage = (c1 + c2) / 2

    def coverage_increment(self):
        self.coverage += 1


class Graph:

    def __init__(self, k):
        self.vertices = {}
        self.k = k

    def add_read(self, read):
        if len(read) < self.k:
            return

        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        for next_kmer_indx in range(1, len(read) - k + 1, 1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx + k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)

            new_edge = Edge(kmer, next_kmer)
            self.vertices[next_kmer].in_edges[kmer] = [new_edge]
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]
            kmer = next_kmer

    def calc_init_edge_coverage(self):

        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(
                    self.vertices[current_vertex].coverage, self.vertices[next_vertex].coverage)

    def visualize(self,graph):
        vis = Digraph(comment='De Brujin genome assemble graph')

        if graph == 'full':
            for k, v in self.vertices.items():

                vis.node(k, label=f'{k}')
                for kk, vv in v.out_edges.items():
                    vis.edge(k, kk, label=f'{vv[0].seq}')
        else:
            for k, v in self.vertices.items():

                vis.node(k, label=f'cov={v.vertex_coverage}')
                for kk, vv in v.out_edges.items():
                    vis.edge(k, kk, label='cov={cov} len={len}'.format(cov=vv[0].edge_coverage, len=len(vv[0].seq)))

        print(vis.source)
        vis.view()
        vis.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='graph_viualization')
    parser.add_argument('-i', help='Name/path to your input fasta file', type=str)
    parser.add_argument('-k', help='k-mer length', default=3, type=int)
    parser.add_argument('-t', help='full/not', default='full', type=str)
    parser.add_argument('-f', help='forward oriented assembly', dest='feature', action='store_true')
    parser.add_argument('-r', help='reverse complement oriented assembly', dest='feature', action='store_false')
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    i, k, t , feature = args.i, args.k, args.t, args.feature
    my_graph = Graph(k)

    with open(i, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if feature == True:
                read = str(record.seq)
            else:
                read = str(record.reverse_complement().seq)
            my_graph.add_read(read)

    my_graph.calc_init_edge_coverage()
    my_graph.visualize(t)