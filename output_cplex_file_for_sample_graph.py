#import MinimumCommonFVS as MCFVS
from MinimumCommonFVS import MinimumCommonFVS as MCFVS
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import random
import os

import matplotlib.pyplot as plt

for_gnm_random_graph = True
for_scale_free = True
for_CACC_net = True

#Sample graph (gnm_random_graph)
if for_gnm_random_graph:
    exch_rate = [5,10,15,20,25]
    tmp = list(map(lambda x:round(x*0.2,2),range(1,4)))

    random.seed(0)
    num_of_nodes = range(100,801,100)
    exch_rate = [5,10,15,20,25]
    num_of_trials = 10

    rand_num = len(num_of_nodes) * num_of_trials * 70000
    random.seed(1)
    rand = random.sample(range(rand_num),rand_num)

    for index, N in enumerate(num_of_nodes):
        M = 2 * N
        list_of_exchange_times = [int(0.01*i*N) for i in exch_rate]
        for l,num_of_times_exchanged in enumerate(list_of_exchange_times):
            print(list_of_exchange_times)
            for trial in range(num_of_trials):
                hoge = True
                while hoge:
                    seed_of_the_1st_graph = random.sample(rand,1)[0]
                    rand.remove(seed_of_the_1st_graph)
                    g1 = nx.gnm_random_graph(
                    N,M,seed=seed_of_the_1st_graph,directed=True
                    )
                    edges_of_g1 = list(g1.edges())
                    g1.remove_edges_from(edges_of_g1)
                    g1.add_edges_from(list(set(edges_of_g1)))
                    seed_of_exchange = random.sample(rand,1)[0]
                    rand.remove(seed_of_exchange)
                    try:
                        g2 = MCFVS.graph_with_2_edges_recombined_N_times(
                        g1,
                        num_of_times_exchanged,
                        seed=seed_of_exchange
                        )
                        graphs = [g1,g2]
                        graph_num = len(graphs)
                        hoge = False
                    except:
                        del seed_of_the_1st_graph
                        del seed_of_exchange
                        hoge = True

                original_edges_of_g1 = g1.number_of_edges()
                original_edges_of_g2 = g2.number_of_edges()
                graphs = MCFVS.graph_editing_for_faster_calculations(graphs)
                g1_of_edges_after_graph_cutting = g1.number_of_edges()
                g2_of_edges_after_graph_cutting = g2.number_of_edges()
                if set(g1.edges) == set(g2.edges):
                    print('In considering FVS, these are the same graphs: seed1;{},seed_of_exchange:{}'.format())
                else:
                    dir_name_for_scale_free_graph = './[random_graph]cplex_files_graph_num[{}]'.format(graph_num)

                    if os.path.isdir(dir_name_for_scale_free_graph+'/node_size[{}]'.format(N)) == False:
                        os.makedirs(dir_name_for_scale_free_graph+'/node_size[{}]'.format(N))

                    file0 = dir_name_for_scale_free_graph + '/node_size[{}]'.format(N) + \
                    '/random_graph_seed_of_first_graph[{}]multi_graphs_nodes[{}]original_edges[{}]exch_rate[{}]exch_num[{}].mod'.\
                    format(seed_of_the_1st_graph,N,original_edges_of_g1,exch_rate[l],num_of_times_exchanged)
                    file1 = dir_name_for_scale_free_graph + '/node_size[{}]'.format(N) + \
                    '/random_graph_seed_of_first_graph[{}]first_graph_nodes[{}]original_edges[{}]edges_after_graph_cutting[{}]exch_rate[{}]exch_num[{}].mod'.\
                    format(seed_of_the_1st_graph,N,original_edges_of_g1,g1_of_edges_after_graph_cutting,exch_rate[l],num_of_times_exchanged)
                    file2 = dir_name_for_scale_free_graph + '/node_size[{}]'.format(N) + \
                    '/random_graph_seed_of_first_graph[{}]exchenged_graph_seed[{}]nodes[{}]original_edges[{}]edges_after_graph_cutting[{}]exch_rate[{}]exch_num[{}].mod'.\
                    format(seed_of_the_1st_graph,seed_of_exchange,N,original_edges_of_g2,g2_of_edges_after_graph_cutting,exch_rate[l],num_of_times_exchanged)
                    filenames = [file0,file1,file2]

                    MCFVS.OutputCPLEXFile(graphs,filenames,Need_MultilayerFVS=True)

##Sample graph (scale_free_graph)
if for_scale_free:
    random.seed(0)
    num_of_nodes = [100,300,500,700,1000,2000,3000,4000,5000,6000]
    exch_rate = [5,10,15,20,25]
    num_of_trials = 10

    list_of_exchange_times = [
    [0.01*i*N for i in exch_rate] for N in num_of_nodes
    ]
    rand_num = len(num_of_nodes) * len(list_of_exchange_times) \
    * num_of_trials * 50000
    random.seed(1)
    rand = random.sample(range(rand_num),rand_num)

    for index, N in enumerate(num_of_nodes):
        list_of_exchange_times = [int(0.01*i*N) for i in exch_rate]
        if N < 1000:
            tmp = list(map(lambda x:round(x*0.2,2),range(1,4)))
            parameters = [(a,b,c) for a in tmp for b in tmp for c in tmp if a+b+c==1]
        elif N >= 1000:
            parameters = [(0.4,0.2,0.4),(0.4,0.4,0.2),('default','default','default')]

        for l,num_of_times_exchanged in enumerate(list_of_exchange_times):
            for abc in parameters:
                for trial in range(num_of_trials):
                    hoge = True
                    while hoge:
                        seed_of_the_1st_graph = random.sample(rand,1)[0]
                        rand.remove(seed_of_the_1st_graph)
                        if abc == ('default','default','default'):
                            g1 = nx.generators.directed.scale_free_graph(
                            N,
                            seed=seed_of_the_1st_graph
                            )
                        else:
                            g1 = nx.generators.directed.scale_free_graph(
                            N,
                            alpha=abc[0],
                            beta=abc[1],
                            gamma=abc[2],
                            seed=seed_of_the_1st_graph
                            )

                        edges_of_g1 = list(g1.edges())
                        g1.remove_edges_from(edges_of_g1)
                        g1.add_edges_from(list(set(edges_of_g1)))
                        seed_of_exchange = random.sample(rand,1)[0]
                        rand.remove(seed_of_exchange)
                        try:
                            g2 = MCFVS.graph_with_2_edges_recombined_N_times(
                            g1,
                            num_of_times_exchanged,
                            seed=seed_of_exchange
                            )
                            graphs = [g1,g2]
                            graph_num = len(graphs)
                            hoge = not MCFVS.scale_free_property_check(graphs)
                        except:
                            del seed_of_the_1st_graph
                            del seed_of_exchange
                            hoge = True

                    original_edges_of_g1 = g1.number_of_edges()
                    original_edges_of_g2 = g2.number_of_edges()
                    graphs = MCFVS.graph_editing_for_faster_calculations(graphs)
                    g1_of_edges_after_graph_cutting = g1.number_of_edges()
                    g2_of_edges_after_graph_cutting = g2.number_of_edges()
                    if set(g1.edges) == set(g2.edges):
                        print('In considering FVS, these are the same graphs: seed1;{},seed_of_exchange:{}'.format())
                    else:
                        dir_name_for_scale_free_graph = './[scale_free]cplex_files_graph_num[{}]'.format(graph_num)
                        if os.path.isdir(dir_name_for_scale_free_graph+'/node_size[{}]'.format(N)) == False:
                            os.makedirs(dir_name_for_scale_free_graph+'/node_size[{}]'.format(N))

                        file0 = dir_name_for_scale_free_graph + '/node_size[{}]'.format(N) + \
                        '/scale_free_abc[{}_{}_{}]seed_of_first_graph[{}]multi_graphs_nodes[{}]original_edges[{}]exch_rate[{}]exch_num[{}].mod'.\
                        format(abc[0],abc[1],abc[2],seed_of_the_1st_graph,N,original_edges_of_g1,exch_rate[l],num_of_times_exchanged)
                        file1 = dir_name_for_scale_free_graph + '/node_size[{}]'.format(N) + \
                        '/scale_free_abc[{}_{}_{}]seed_of_first_graph[{}]first_graph_nodes[{}]original_edges[{}]edges_after_graph_cutting[{}]exch_rate[{}]exch_num[{}].mod'.\
                        format(abc[0],abc[1],abc[2],seed_of_the_1st_graph,N,original_edges_of_g1,g1_of_edges_after_graph_cutting,exch_rate[l],num_of_times_exchanged)
                        file2 = dir_name_for_scale_free_graph + '/node_size[{}]'.format(N) + \
                        '/scale_free_abc[{}_{}_{}]seed_of_first_graph[{}]exchenged_graph_seed[{}]nodes[{}]original_edges[{}]edges_after_graph_cutting[{}]exch_rate[{}]exch_num[{}].mod'.\
                        format(abc[0],abc[1],abc[2],seed_of_the_1st_graph,seed_of_exchange,N,original_edges_of_g2,g2_of_edges_after_graph_cutting,exch_rate[l],num_of_times_exchanged)
                        filename = [file0,file1,file2]
                        MCFVS.OutputCPLEXFile(graphs,filename,Need_MultilayerFVS=True)

#CACC network
if for_CACC_net:
    DirNameforRealBioData = \
    './Correspondence_table_of_CACC_net_for_cplex'
    if os.path.isdir(DirNameforRealBioData) == False:
        os.makedirs(DirNameforRealBioData)
    MCFVS._Generate_ids_for_cplex_of_graph_nodes()
