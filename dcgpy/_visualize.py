def _dcgpann_visualize(self, show_connections = True, fill_color = 'w', show_nonlinearities = False, active_connection_alpha = 0.1, inactive_connection_alpha = 0.01, legend = True, axes = None):
    """visualize(show_connections = True, fill_color = 'w', show_nonlinearities = False, active_connection_alpha = 0.1, inactive_connection_alpha = 0.01, legend = True)

    Visualizes the dCGPANN expression

    Args:
        show_connections (``bool``): shows all the connections between nodes (for large dCGPANN can be quite time consuming and cluttered)
        fill_color (``bool``): the fill color of all nodes
        show_nonlinearities (``bool``): shows also what nonlinearity is where and adds a legend
        active_connection_alpha (``bool``): the alpha used to show active connections
        inactive_connection_alpha (``bool``): the alpha used to show inactive connections
        legend (``bool``): shows a legend when show_nonlinearities is also True



    Examples:

    >>> from dcgpy import *
    >>> nonlinearities = dcgpy.kernel_set_double(["sig", "ReLu", "tanh"])
    >>> dcgpann = dcgpy.expression_ann_double(inputs=3, outputs=4, rows=15, cols=8, 
                                      levels_back=2, arity=4, kernels=nonlinearities(), seed=32)
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
            pos[n+r*c+inp] = [xpos, ypos + dr*inp + offset]      

    # ------------------------------_HERE WE PLOT -----------------------------------------------
    if show_connections:
        # We loop over all cartesian nodes and plot intranodes connections
        weight_normalization = max(np.abs(self.get_weights()))
        for node_id in range(n, n+r*c):
            start = gene[node_id]
            connections = x[start+1:start+1+self.get_arity(node_id)]
            for j, target_node in enumerate(connections):
                if node_id in an:
                    alpha = abs(self.get_weight(node_id, j)) / weight_normalization * active_connection_alpha
                else:
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