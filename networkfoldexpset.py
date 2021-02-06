#!/usr/bin/env python

from base64 import b64encode
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import time
import math
from granatum_sdk import Granatum
import networkx as nx
import re
from networkx.drawing.nx_agraph import write_dot
import os

def geturl(gene_id):
    return "https://www.genecards.org/cgi-bin/carddisp.pl?gene={}".format(gene_id)

def main():
    tic = time.perf_counter()

    gn = Granatum()

    clustersvsgenes = gn.pandas_from_assay(gn.get_import('clustersvsgenes'))
    max_dist = gn.get_arg('max_dist')
    min_zscore = gn.get_arg('min_zscore')

    clustercomparisonstotest = list(clustersvsgenes.index)

    G = nx.MultiDiGraph()
    clusternames = list(clustersvsgenes.T.columns)
    individualclusters = [n[:n.index(" vs rest")] for n in clusternames if n.endswith("vs rest")]
    print(individualclusters, flush=True)
    for cl in individualclusters:
        G.add_node(cl)

    # {pathway : {"cluster1":score1, "cluster2":score2}, pathway2 : {}}
    resultsmap = {}
    relabels = {}
    keys = {}
    currentkeyindex = 0
    for gene_id in clustersvsgenes.columns: 
        for cluster in clustercomparisonstotest:
            score = clustersvsgenes.loc[cluster, gene_id]
            if score >= min_zscore:
                if not gene_id in keys:
                    # First check if within distance of another group 
                    closestkey = None
                    closestkeyvalue = 1.0e12
                    for key in keys:
                        gene_values = clustersvsgenes.loc[:, gene_id]
                        ref_values = clustersvsgenes.loc[:, key]
                        sc = np.sqrt(np.nansum(np.square(gene_values-ref_values)))
                        if sc <= max_dist and sc < closestkeyvalue:
                            closestkeyvalue = sc
                            closestkey = key
                    if closestkey == None:
                        keys[gene_id] = currentkeyindex + 1
                    else:
                        keys[gene_id] = keys[closestkey]
                            
                print("Score = {}".format(score), flush=True)
                olddict = resultsmap.get(gene_id, {})
                olddict[cluster] = score
                resultsmap[gene_id] = olddict
                from_to = re.split(' vs ', cluster)
                if from_to[1] != 'rest':
                    G.add_weighted_edges_from([(from_to[1], from_to[0], score*2.0)], label=str(keys[gene_id]), penwidth=str(score*2.0))
                else:
                    relabel_dict = relabels.get(from_to[0], "")
                    if relabel_dict == "":
                        relabel_dict = from_to[0] + ": " + str(keys[gene_id])
                    else:
                        relabel_dict = relabel_dict + ", " + str(keys[gene_id])
                    relabels[from_to[0]] = relabel_dict
                currentkeyindex = max(currentkeyindex, keys[gene_id])

    print("Relabels {}".format(relabels), flush=True)
    G = nx.relabel_nodes(G, relabels)
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G, 'label')
    write_dot(G, 'plot.dot')
    os.system("dot plot.dot -Tpng -Gdpi=600 > plot.png")
    with open('plot.png', "rb") as f:
        image_b64 = b64encode(f.read()).decode("utf-8")

    gn.results.append(
        {
                "type": "png",
                "width": 650,
                "height": 480,
                "description": 'Network of clusters based on expression',
                "data": image_b64,
            })

    footnote = ""
    inv_map = {}
    for k, v in keys.items():
        inv_map[v] = inv_map.get(v, []) + [k]

    for k, v in sorted(inv_map.items(), key=lambda item: item[0]):
        newv = map(lambda gene: "[{}]({})".format(gene, geturl(gene)), v)
        vliststr = ", ".join(newv)
        newstr = "{}: {}".format(v, vliststr)
        if footnote == "" :
            footnote = newstr
        else:
            footnote = footnote + "  \n"+newstr

    gn.add_result(footnote, "markdown")

    # gn.export(return_df.T.to_csv(), 'differential_gene_sets.csv', kind='raw', meta=None, raw=True)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
