def _dcgpann_visualize(self, show_connections = True, fill_color = 'w', active_connection_alpha = 0.1, inactive_connection_alpha = 0.01, legend = True, axes = None, show_inactive = False, show_nonlinearities = False):
    """visualize(show_connections = True, fill_color = 'w', show_nonlinearities = False, active_connection_alpha = 0.1, inactive_connection_alpha = 0.01, legend = True)

    Visualizes the dCGPANN expression

    Args:
        show_connections (``bool``): shows active connections between nodes
        show_inactive (```bool```): shows also inactive connections between nodes
        active_connection_alpha (``bool``): the alpha used to show active connections
        inactive_connection_alpha (``bool``): the alpha used to show inactive connections
        fill_color (``str`` or ```RGB values```): the fill color of all nodes
        show_nonlinearities (``bool``): color codes the nodes with the contained nonlinearity
        legend (``bool``): shows a legend (only when show_nonlinearities is also True)



    Examples:

    >>> from dcgpy import *
    >>> nonlinearities = dcgpy.kernel_set_double(["sig", "ReLu", "tanh"])
    >>> dcgpann = dcgpy.expression_ann_double(inputs=3, outputs=4, rows=15, cols=8, levels_back=2, arity=4, kernels=nonlinearities(), seed=32)
    >>> dcgpann.randomise_weights()
    >>> dcgpann.visualize(show_nonlinearities=False)
    """



    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Failed to import the required module matplotlib")
        raise
        
    try:
        import numpy as np
    except ImportError:
        print("Failed to import the required module numpy")
        raise

    # We prepare the plot
    if axes is None:
        fig, ax = plt.subplots() 
    else:
        ax = axes
    plt.axis('off')
    ax.grid(False)
    color_map = ['k', 'b','g','r','c','m', 'y']

    # We get some properties of the dCGPANN expression
    n = self.get_n()
    m = self.get_m()
    r = self.get_rows()
    c = self.get_cols()
    gene = self.get_gene_idx()

    # We get the active nodes id and add the output nodes id (they do not have an id in CGP encoding)
    an = self.get_active_nodes()
    an = an + [it for it in range(n+r*c, n+r*c+m)]

    # We get the cromosome
    x = self.get()

    # Grid parameters: dr is the separation between rows and dc between columns. radius is the neuron radius
    dr = 1./(r-1)
    dc = 1./(c-1)
    radius = 1./max(r,c)/4.

    # pos will hold the nodes poistions
    pos = {}
    nodeid = n
    # Cartesian nodes positions
    for col in range(c):
        xpos = 0.
        for row in range(r):
            ypos = 1.
            pos[nodeid] = [xpos + dc*col, ypos - dr*row]
            nodeid+=1
    # Input nodes positions
    xpos = - dc 
    ypos = 0
    offset = (1. - (n-1) * dr) / 2.
    for inp in range(n):
            pos[inp] = [xpos, ypos + dr*inp + offset]
    # Output nodes positions
    xpos = 1 + dc 
    offset = (1. - (m-1) * dr) / 2.
    for inp in range(m):
            pos[n+r*c+inp] = [xpos, ypos + (dr*(m-1)) + offset - dr*inp]      

    # ------------------------------_HERE WE PLOT -----------------------------------------------
    if show_connections:
        # We loop over all cartesian nodes and plot intranodes connections
        weight_normalization = max(np.abs(self.get_weights()))
        for node_id in range(n, n+r*c):
            start = gene[node_id]
            connections = x[start+1:start+1+self.get_arity(node_id)]
            for j, target_node in enumerate(connections):
                if node_id in an:
                    alpha = abs(self.get_weight(node_id, j)) / weight_normalization * active_connection_alpha + 0.1
                    ax.plot([pos[node_id][0], pos[target_node][0]],[pos[node_id][1], pos[target_node][1]],alpha=alpha, color='k', zorder=-10)
                elif show_inactive:
                    alpha = inactive_connection_alpha
                    ax.plot([pos[node_id][0], pos[target_node][0]],[pos[node_id][1], pos[target_node][1]],alpha=alpha, color='k', zorder=-10)
                
        # And the connections to the output
        for i in range(m):
            start = n+r*c
            ax.plot([pos[start+i][0], pos[x[-m+i]][0]],[pos[start+i][1], pos[x[-m+i]][1]],alpha=active_connection_alpha, color='k', zorder=-10)

    # We plot the cartesian nodes
    for node_id in pos:
        if node_id < n and node_id >= n+r*c and show_nonlinearities:
            # For input and output nodes we use gray as a filler when visualizing nonlinearities
            circle = plt.Circle((pos[node_id][0], pos[node_id][1]), radius, color='gray', fill = True)
        else:
            circle = plt.Circle((pos[node_id][0], pos[node_id][1]), radius, color=fill_color, fill = True)
        if node_id in an:
            alpha = 0.8
        else:
            alpha = 0.1
        ax.add_artist(circle)
        if show_nonlinearities and node_id >=n and node_id < n+r*c:
            nl_type = x[gene[node_id]]
            nl_color = color_map[nl_type % len(color_map)]
        else:
            nl_color = 'k'
        circle = plt.Circle((pos[node_id][0], pos[node_id][1]), radius, color=nl_color, alpha=alpha,fill = False)
        ax.add_artist(circle)

    xmin = min([pos[key][0] for key in pos]) - radius * 2
    xmax = max([pos[key][0] for key in pos]) + radius * 2
    ymin = min([pos[key][1] for key in pos]) - radius * 2
    ymax = max([pos[key][1] for key in pos]) + radius * 2

    # We plot a legend if we are visualising nonlinearities
    if show_nonlinearities and legend:
        legend_x = xmax + 0.3
        legend_y = ymax - radius * 3
        sep = max(3 *radius, 0.08)
        for i,kernel in enumerate(self.get_f()):
            circle = plt.Circle((legend_x, legend_y - i * sep), radius, color=color_map[i % len(color_map)], alpha=0.8,fill = False)
            plt.text(legend_x + 3 * radius, legend_y - i *sep, kernel.__repr__())
            ax.add_artist(circle)
        xmax = xmax + 0.4

    # We set the axis
    ax.set_ylim((ymin, ymax))
    ax.set_xlim((xmin, xmax))

    return ax

def _graphviz_visualize(self, in_sym = [], draw_inactive = True, draw_weights = False):
    """
    visualize(in_sym = [], draw_inactive = True, draw_weights = False)
    Visualizes the d-CGP expression
    Visualizes the graph of the d-CGP expression
    Note:
        This method requires the `graphviz`` module installed in your Python system
    Args:
        in_sym (a ``List[str]``): input symbols. Its length must either match the number of inputs or be zero (to visualize them as x_i)
        draw_inactive (a ``bool``): indicates whether to draw inactive nodes
        draw_weights (a ``bool``): indicates whether to draw connection weights symbols
    Returns:
        The ``graphviz.Digraph`` for the given expression
    Raises:
        ImportError: if module pygraphviz is not installed in your Python system
        ValueError: if in_sym is nonempty but its length does not match the number of inputs
    Examples:
        >>> ex = dcgpy.expression_double(2,1,3,3,2,2,dcgpy.kernel_set_double(["sum","diff"])(),0)
        >>> graph = ex.visualize(['x', 'c'], True, False)
    """

    from graphviz import Digraph

    n = self.get_n()

    if len(in_sym) != 0 and len(in_sym) != n:
        raise ValueError("The length of in_sym must either match the number of inputs or be zero")

    x = self.get()
    m = self.get_m()
    r = self.get_rows()
    c = self.get_cols()
    f = self.get_f()
    active_nodes = self.get_active_nodes()

    # bool vector of active nodes
    is_active = [False] * (n + r * c)
    for i in range(len(active_nodes)):
        is_active[active_nodes[i]] = True

    G = Digraph(graph_attr={'rankdir': 'LR'})
    
    # force the nodes to be placed in the right ranks
    inputs = Digraph(graph_attr={'rank': 'same'})
    for i in range(n):
        inputs.node('n' + str(i))
    G.subgraph(inputs)
    for i in range(c):
        col = Digraph(graph_attr={'rank': 'same'})
        for j in range(r):
            col.node('n' + str(n + (i * r) + j))
        G.subgraph(col)
        if i == 0:
            for j in range(n):
                for k in range(r):
                    G.edge('n' + str(j),'n' + str(n + k), style = 'invis')
        else:
            for j in range(r):
                for k in range(r):
                    G.edge('n' + str(n + (i - 1) * r + j),'n' + str(n + i * r + k), style = 'invis')
    outputs = Digraph(graph_attr={'rank': 'same'})
    for i in range(m):
        outputs.node('n' + str(n + r * c + i))
    G.subgraph(outputs)
    for i in range(r):
        for j in range(m):
            G.edge('n' + str(n + (c - 1) * r + i),'n' + str(n + r * c + j), style = 'invis')

    # input nodes
    for i in range(n):
        if len(in_sym) != 0:
            xlabel = in_sym[i]
        else:
            xlabel = '<x<sub>' + str(i) + '</sub>>'
        G.node('n' + str(i), label = xlabel, shape = 'circle', style = 'bold')

    # function nodes and connections
    for i in range(r * c):
        node_id = n+i
        if is_active[node_id]:
            nstyle = 'solid'
            estyle = 'solid'
            col = 'black'
        elif draw_inactive:
            nstyle = 'dashed'
            estyle = 'dotted'
            col = 'grey70'
        else:
            nstyle = 'invis'
            estyle = 'invis'
            col = 'black'
        op = str(f[x[self.get_gene_idx()[node_id]]])
        if op == 'sum':
            op = '+'
        elif op == 'diff':
            op = '-'
        elif op == 'mul':
            op = '*'
        elif op == 'div' or op == 'pdiv':
            op = '/'
        G.node('n' + str(node_id), label =  op, shape = 'circle', style = nstyle, color = col, fontcolor = col)
        for j in range(self.get_arity(node_id)):
            if j == 0:
                ah = 'lnormal'
            elif j == 1:
                ah = 'rnormal'
            else:
                ah = 'normal'
            if draw_weights:
                elabel = '<w<sub>' + str(node_id) + ',' + str(j) + '</sub>>'
            else:
                elabel = ''
            G.edge('n' + str(x[self.get_gene_idx()[node_id] + 1 + j]), 'n' + str(node_id), label = elabel, arrowhead = ah, style = estyle, color = col, fontcolor = col)

    # output nodes
    for i in range(m):
        G.node('n' + str(n + r * c + i), label = '<o<sub>' + str(i) + '</sub>>', shape = 'circle', style = 'bold')
        G.edge('n' + str(x[len(x) - m + i]), 'n' + str(n + r * c + i))

    return G
