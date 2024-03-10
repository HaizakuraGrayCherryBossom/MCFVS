from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random
import time
import math
import csv
import re
import os

#graphs = [g1,g2,...]. then the program runs based on the nodes in "g1".
class MinimumCommonFVS():
    #Output edge list created by randomly replacing edges
    #g: a graph (of networkx), Number_of_times_to_replace: number of times to replace
    def graph_with_2_edges_recombined_N_times(g,Number_of_times_to_replace,seed=None):
        old, new = [], []
        edge_list = list(g.edges())
        original_edge_list, edge_list = [edge_list] * 2
        i, i_max = 1, int(math.factorial(len(edge_list))/(math.factorial(len(edge_list)-2)*2))
        while len(new) < 2 * Number_of_times_to_replace:
            if i <= i_max:
                random.seed(seed)
                edge1, edge2 = random.sample(edge_list,2)
                random.shuffle(edge_list)
                tuple1, tuple2 = (edge1[0],edge2[1]), (edge2[0],edge1[1])
                if not (tuple1 in edge_list+old+new) and not(tuple2 in edge_list+old+new):
                    edge_list.remove(edge1), edge_list.remove(edge2), \
                    old.append(edge1), old.append(edge2), \
                    new.append(tuple1), new.append(tuple2)
            else:
                raise ValueError('\033[31m'+'There is no possible branch exchange.'+'\033[0m')
            i += 1
        new_edge_list = [i for i in original_edge_list if i not in old] + new
        g2 = nx.DiGraph()
        g2.add_nodes_from(list(g.nodes))
        g2.add_edges_from(new_edge_list)
        if len(g.edges)!=len(g2.edges):
            raise ValueError('The number of edges before and after replacement does not match. Are there any problems such as overlapping edges in the input graph?')
        return g2

#i:int >= 2 (Number of edges in a cycle
#(each cycle is detected as 2,3,4 and so on up to the maximum value))
    def _detect_cycles_with_time_limit(nodes_with_self_loop,G):
        tmp = [nodes_with_self_loop]
        n_gon = 1
        while n_gon == 1 or n_gon != 13:
            n_gon += 1
            start = time.time()
            exec('all_the_nodes = G[list(G.keys())[0]];'\
            +'ll = G.values();tmp.append(list(set([tuple(sorted(['+\
            ','.join(['node{}'.format(j) for j in range(1,n_gon+1)])+'])) '\
            +'for graph in ll for node1 in all_the_nodes '\
            + ''.join(['for node{} in [node{} for node{} '.format(k,k,k)\
            +'in set(graph.successors(node{})) if '.format(k-1)+' and '\
            .join(['node{}!=node{}'.format(k,j) for j in range(1,k)])+'] ' \
            if k!=n_gon else 'for node{} in [node{} for node{} '.format(k,k,k)\
            +'in set(graph.successors(node{})) if '.format(k-1)+' and '\
            .join(['node{}!=node{}'.format(k,j) for j in range(1,k)]) \
            for k in range(2,n_gon+1)])\
            +' and (node1 in set(graph.successors(node{})))] if '.format(n_gon)\
            +'all([1 if not set(n_tuple)<set(['\
            +','.join(['node{}'.format(j) for j in range(1,n_gon+1)])\
            +']) else 0 for n_tuple in {}])'.format(sum(tmp,[]))\
            +'])))',{'G':G,'tmp':tmp}\
            )
            end = time.time()
            if end - start > 5:
                break
        return [\
        ' + '.join(['x{} '.format(k) for k in j]) + '>= 1\n' \
        for n_gon in tmp[1:] for j in n_gon\
        ]

    #G: graphs (dict) <ex> {id1:g1,id2:g2,...}, filename: file name (str)
    #Perform ILP formulation from GRAPH data. Then CPLEX file output.
    def _ILP_formulation_and_CPLEX_file_output(G,filename,constraints_cycle_detection):
        graph_ids = list(G.keys())
        nodes = G[list(G.keys())[0]]
        N = len(nodes)
        display = 'display solution variables -'
        obj = ''.join([' x{} +'.format(i) for i in nodes])[1:-2]
        x = obj.replace('+ ','')
        sbj_list = []

        #Subject to: Self loop
        self_loop_edge_candidates = {(i,i) for i in nodes}
        nodes_with_self_loop = list(\
        {(edge[0],) for graph in G.values() \
        for edge in self_loop_edge_candidates \
        if edge in graph.edges()}\
        )
        sbj_list += [\
        'x{0} <= 1\n'.format(sl[0])+'x{0} >= 1\n'.format(sl[0]) \
        for sl in nodes_with_self_loop\
        ]
        
        #Subject to: Cycle detection
        if constraints_cycle_detection:
            tmp = MinimumCommonFVS._detect_cycles_with_time_limit(nodes_with_self_loop,G)
            sbj_list += tmp

        #Subject to: No destinations containing isolated nodes.
        no_destination = [all([set(G[v].successors(node)) == set() \
        for v in graph_ids if v in G]) for node in nodes]
        sbj_list += ['x{} >= 0\n'.format(node)+'x{} <= 0\n'.format(node) \
        for index,node in enumerate(nodes) if no_destination[index]]
        #Subject to: No the past containing isolated nodes.
        no_the_past = [all([set(G[v].predecessors(node)) == set() \
        for v in graph_ids if v in G]) for node in nodes]
        sbj_list += ['x{} >= 0\n'.format(node)+'x{} <= 0\n'.format(node) \
        for index,node in enumerate(nodes) if no_the_past[index] \
        ]
        #Subject to:conditional formula for k_parameters
        sbj_list += sum(\
        [['k({0},{1}) - k({0},{2}) + {3} x{1} >= 1\n'.format(h,edge[0],edge[1],N) \
        for edge in G[v].edges() if edge[0]!=edge[1]] for h,v in enumerate(graph_ids) if v in G\
        ],[]\
        )
        #Stringing conditionals
        sbj = ''.join(list(set(sbj_list)))
        #k parameters: CPLEX;generary
        k_parameters = [['k({},{}) '.format(h,node),'0 <= k({},{}) <= {}\n'.format(h,node,N-1)] \
        for h, v in enumerate(graph_ids) for node in nodes.nodes() if (v in G) and ('k({},{})'.format(h,node) in sbj)]
        k, k_b = ''.join([k[0] for k in k_parameters]),''.join([k[1] for k in k_parameters])

        #Write ILP formula to file and output
        with open(filename,'w') as f:
            f.write('enter {}'.format(filename.split('/')[-1])[:-4]+'\n')#[:-4]:Remove the extension
            f.write("minimize\n")
            f.write(obj+'\n')
            f.write('subject to\n')
            f.write(sbj)
            f.write('bounds\n')
            f.write(k_b)
            f.write('generary\n')
            f.write(k+'\n')
            f.write('binary\n')
            f.write(x+'\n')
            f.write('end\n')
            f.write('optimize\n')
            f.write(display+'\n')
            f.write('quit')

#Output CPLEX file:
#G is a list of graphs (list), and "filename" is, for example, ['dirname/multi_hoge.mod','dirname/hoge1.mod ','dirname/hoge2.mod',...] (list)
#Filenames must be as follows: [graph of multilayer network, filename of graph1, filename of graph2, ...]'
    def OutputCPLEXFile(G,filename,constraints_cycle_detection=True,Need_MultilayerFVS=None):
        #If the filename is a single str type, make it a list with one element.
        if type(filename) != list:
            if type(filename) == str:
                filename = [filename]
            else:
                raise ValueError('The filename when considering only multi-layer networks must be as follows: "graph of multilayer network" or ["graph of multilayer network"]')

        #Preparation for multilayer networks
        graphs = {i:g for i, g in enumerate(G)}
        MinimumCommonFVS._ILP_formulation_and_CPLEX_file_output(graphs,filename[0],constraints_cycle_detection)

        if Need_MultilayerFVS:
            #Error when filename and number of graphs do not match.
            if len(filename) != len(G) + 1:
                raise ValueError('Filenames must be as follows: [graph of multilayer network, filename of g1, filename of g2, ...]')
            for i, g in enumerate(G):
                MinimumCommonFVS._ILP_formulation_and_CPLEX_file_output({i:g},filename[i+1],constraints_cycle_detection)

#graphs:nx.Graph()
    def graph_editing_for_faster_calculations(graphs):
        #N = len(graphs[list(graphs.keys())[0]])
        N = graphs[0].number_of_nodes()
        self_loop_edge_candidates = {(i,i) for i in range(N)}
        self_loops = list(set([\
        edge for graph in graphs \
        for edge in self_loop_edge_candidates if edge in graph.edges() \
        ]))
        [\
        (graph.remove_edges_from(list(graph.in_edges(edge[0]))),\
        (graph.remove_edges_from(list(graph.edges(edge[0]))))) \
        for graph in graphs for edge in self_loops\
        ], [graph.add_edges_from(self_loops) for graph in graphs]
        dummy = 0
        edges_of_graph0, edges_of_graph0_new = [], []
        while edges_of_graph0 != edges_of_graph0_new or dummy == 0:
            edges_of_graph0 = set(graphs[0].edges())
            [\
            ([graph.remove_edges_from(list(graph.in_edges(node))) for graph in graphs \
            if all([(set(graph.successors(node)) == set()) for graph in graphs])\
            ], [\
            graph.remove_edges_from(list(graph.edges(node))) for graph in graphs \
            if all([set(graph.predecessors(node)) == set() for graph in graphs])\
            ]) for node in graphs[0].nodes()\
            ]
            edges_of_graph0_new = set(graphs[0].edges())
            dummy += 1
        return graphs

    def scale_free_property_check(graphs):
        if type(graphs) != list:
            graphs = [graphs]
        f = lambda x, c, l: l * (x**(-c))
        check = []
        for g in graphs:
            degree_ = [d for i,d in g.degree()]
            k = sorted(list(set(degree_)))
            p_k = [degree_.count(i)/len(degree_) for i in k]
            popt, pcov = curve_fit(f,k,p_k)
            check.append(1.5 < popt[0] < 4.5)
        if all(check):
            return True
        else:
            return False
