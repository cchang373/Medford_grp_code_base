#!/usr/bin/env python276
import cPickle as pickle

# Import graphviz
import gv

# Import pygraph
from pygraph.classes.graph import graph
from pygraph.classes.digraph import digraph
from pygraph.algorithms.searching import breadth_first_search, depth_first_search
from pygraph.algorithms.minmax import shortest_path
from pygraph.readwrite.dot import write

def make_graph(rxns,name='',subscript_characters = ['x','y','z','w','2','3','4'],node_arg_dict={}):
    nodes = []
    edges = []
    edge_args = {}
    default_edge_args = {'penwidth':'3'}
    for i,rxn in enumerate(rxns):
        if len(rxn) > 2:
            rxn_args = rxn[-1]
        else:
            rxn_args = {}

        nodes.append(str(i))
        node_arg_dict[str(i)] = rxn_args

        if '' in rxn[0]:
            rxn[0].remove('')
        for ri in rxn[0]:
            if ri not in nodes and ri:
                nodes.append(ri)
            if [ri,str(i)] not in edges:
                edges.append([ri,str(i)])
                edge_args[str([ri,str(i)])] = rxn_args

        if '' in rxn[1]:
            rxn[1].remove('')
        for ri in rxn[1]:
            if ri not in nodes and ri:
                nodes.append(ri)
            if [str(i),ri] not in edges:
                edges.append([str(i),ri])
                edge_args[str([str(i),ri])] = rxn_args


    # Graph creation
    default_node_args = {'color':'gray', 'bgcolor':'0000ff80', 'penwidth':'2'}
    gr = digraph()
    for n in nodes:
        node_args = default_node_args.copy()
        node_args.update(node_arg_dict.get(n,{}))

        try:
            int(n)
            node_args['label'] = '.'
            node_args['fixedsize'] = 'true'
            node_args['width'] = 0.1
            node_args['height'] = 0.1
            attrs = node_args.items()
            gr.add_node(n,attrs=attrs)
        except:
            nlab = n
            for ch in subscript_characters:
                if ch in n:
                    nlab = nlab.replace(ch, '<sub>'+ch+'</sub>')
            if nlab != n:
                node_args['label'] = '<'+nlab+'>'

            if '(g)' in n:
                node_args['shape'] = 'invhouse'

                if 'CH3CH2OH' in n:
                    node_args['color'] = 'red'
                elif 'CH3OH' in n:
                    node_args['color'] = 'green'
                elif 'CH4' in n:
                    node_args['color'] = 'blue'
                else:
                    node_args['color'] = 'black'

                node_args['bgcolor'] = '00ff0080'
                attrs = node_args.items()
                gr.add_node(n,attrs=attrs)
            else:
                attrs = node_args.items()
                gr.add_node(n,attrs=attrs)

    for edge in edges:
        edge_i_args = default_edge_args.copy()
        if str(edge) in edge_args:
            edge_i_args.update(edge_args[str(edge)])
        attrs = edge_i_args.items()
        gr.add_edge(edge,attrs=attrs)

    dot = write(gr)
    gvv = gv.readstring(dot)
    gv.setv(gvv,'name',name)
    gv.setv(gvv,'size','8,11')
    gv.layout(gvv,'dot')
    gv.render(gvv,'pdf',name+'_network.pdf')


if __name__ == '__main__':
    rxns = [
            [['CO(g)'],['CO*'],{'color':'black'}],
            [['H2(g)'],['H*'],{'color':'black'}],
            [['CO*','H*'],['CHxOHy*'],{'color':'gray'}],
            [['CHxOHy*','H*'],['CHxOHy*'],{'color':'gray'}],
            [['CHxOHy*','CHxOHy*'],['HxOCHyCHzOHw*'],{'color':'cyan'}],
            [['H*','HxOCHyCHzOHw*'],['HxOCHyCHzOHw*'],{'color':'gray'}],
            [['CHxOHy*'],['CHx*','OHx*'],{'color':'orange'}],
            [['H*','CHx*'],['CHx*'],{'color':'gray'}],
            [['CHxOHy*','CHx*'],['CHxCHyOHz*'],{'color':'cyan'}],
            [['H*','CHxCHyOHz*'],['CHxCHyOHz*'],{'color':'gray'}],
#            [['CHx*','CHx*'],['CHxCHy*'],{'color':'cyan'}],
            [['HxOCHyCHzOHw*'],['CHxCHyOHz*','OHx*'],{'color':'orange'}],
#            [['CHxCHyOHz*'],['CHxCHy*','OHx*'],{'color':'orange'}],
            [['OHx*'],['H2O(g)'],{'color':'black'}],
            [['H*','OHx*'],['OHx*'],{'color':'gray'}],
            [['CHx*'],['CH4(g)'],{'color':'black'}],
            [['CHxOHy*'],['CH3OH(g)'],{'color':'black'}],
#            [['CHxCHy*'],['CH3CH3(g)'],{'color':'black'}],
            [['CHxCHyOHz*'],['CH3CH2OH(g)'],{'color':'black'}],
#            [['CHxCHyOHz*'],['CH3CHO(g)'],{'color':'black'}],
#            [['CO*','OHx*'],['CO2(g)'],{'color':'black'}],
            ]

    make_graph(rxns,'EtOH_OCCO')
