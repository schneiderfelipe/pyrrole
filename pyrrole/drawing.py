#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""Tools for drawing chemical diagrams."""

# TODO: revisit all examples in this module.

import os
import collections as _collections

import numpy as _np
import networkx as _nx
import matplotlib as _mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    _mpl.use('Agg')
import matplotlib.cbook as _cb  # noqa
import matplotlib.pyplot as _plt  # noqa
from matplotlib.colors import Colormap as _Colormap  # noqa
from matplotlib.colors import colorConverter as _colorConverter  # noqa
from matplotlib.collections import LineCollection as _LineCollection  # noqa


def _longest_common_subsequence(x, y):
    """
    Return the longest common subsequence between two sequences.

    Parameters
    ----------
    x, y : sequence

    Returns
    -------
    sequence
        Longest common subsequence of x and y.

    Examples
    --------
    >>> _longest_common_subsequence("AGGTAB", "GXTXAYB")
    ['G', 'T', 'A', 'B']
    >>> _longest_common_subsequence(["A", "GA", "G", "T", "A", "B"],
    ...                             ["GA", "X", "T", "X", "A", "Y", "B"])
    ['GA', 'T', 'A', 'B']

    """
    m = len(x)
    n = len(y)

    # L[i, j] will contain the length of the longest common subsequence of
    # x[0..i - 1] and y[0..j - 1].
    L = _np.zeros((m + 1, n + 1), dtype=int)

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                continue
            elif x[i - 1] == y[j - 1]:
                L[i, j] = L[i - 1, j - 1] + 1
            else:
                L[i, j] = max(L[i - 1, j], L[i, j - 1])

    ret = []

    i, j = m, n
    while i > 0 and j > 0:
        # If current character in x and y are same, then current character is
        # part of the longest common subsequence.
        if x[i - 1] == y[j - 1]:
            ret.append(x[i - 1])
            i, j = i - 1, j - 1
        # If not same, then find the larger of two and go in the direction of
        # larger value.
        elif L[i - 1, j] > L[i, j - 1]:
            i -= 1
        else:
            j -= 1

    return ret[::-1]


def tower_layout(graph, height='freeenergy', scale=None, center=None, dim=2):
    """
    Position all nodes of graph stacked on top of each other.

    Parameters
    ----------
    graph : `networkx.Graph` or `list` of nodes
        A position will be assigned to every node in graph.
    height : `str` or `None`, optional
        The node attribute that holds the numerical value used for the node
        height. This defaults to ``'freeenergy'``. If `None`, all node heights
        are set to zero.
    scale : number, optional
        Scale factor for positions.
    center : array-like, optional
        Coordinate pair around which to center the layout. Default is the
        origin.
    dim : `int`
        Dimension of layout. If `dim` > 2, the remaining dimensions are set to
        zero in the returned positions.

    Returns
    -------
    pos : mapping
        A mapping of positions keyed by node.

    Examples
    --------
    >>> from pyrrole import ChemicalSystem
    >>> from pyrrole.atoms import create_data, read_cclib
    >>> from pyrrole.drawing import tower_layout
    >>> data = create_data([
    ...     read_cclib("data/acetic_acid.out", "AcOH(g)"),
    ...     read_cclib("data/acetic_acid@water.out", "AcOH(aq)")])
    >>> digraph = (ChemicalSystem("AcOH(g) <=> AcOH(aq)", data)
    ...            .to_digraph())
    >>> layout = tower_layout(digraph)
    >>> layout['AcOH(g)']
    array([   0.        , -228.56450866])

    Passing ``scale=1`` means scaling positions to ``(-1, 1)`` in all axes:

    >>> layout = tower_layout(digraph, scale=1)
    >>> layout['AcOH(g)'][1] <= 1.
    True

    """
    # TODO: private function of packages should not be used.
    graph, center = _nx.drawing.layout._process_params(graph, center, dim)

    num_nodes = len(graph)
    if num_nodes == 0:
        return {}
    elif num_nodes == 1:
        return {_nx.utils.arbitrary_element(graph): center}

    paddims = max(0, (dim - 2))

    if height is None:
        y = _np.zeros(len(graph))
    else:
        y = _np.array([data for node, data in graph.nodes(data=height)])
    pos_arr = _np.column_stack([_np.zeros((num_nodes, 1)), y,
                                _np.zeros((num_nodes, paddims))])

    if scale is not None:
        pos_arr = _nx.drawing.layout.rescale_layout(pos_arr,
                                                    scale=scale) + center
    pos = dict(zip(graph, pos_arr))

    # TODO: make test
    return pos


def diagram_layout(graph, height='freeenergy', sources=None, targets=None,
                   pos=None, scale=None, center=None, dim=2):
    """
    Position nodes such that paths are highlighted, from left to right.

    Parameters
    ----------
    graph : `networkx.Graph` or `list` of nodes
        A position will be assigned to every node in graph.
    height : `str` or `None`, optional
        The node attribute that holds the numerical value used for the node
        height. This defaults to ``'freeenergy'``. If `None`, all node heights
        are set to zero.
    sources : `list` of `str`
        All simple paths starting at members of `sources` are considered.
        Defaults to all nodes of graph.
    targets : `list` of `str`
        All simple paths ending at members of `targets` are considered.
        Defaults to all nodes of graph.
    pos : mapping, optional
        Initial positions for nodes as a mapping with node as keys and
        values as a coordinate `list` or `tuple`. If not specified (default),
        initial positions are computed with `tower_layout`.
    scale : number, optional
        Scale factor for positions.
    center : array-like, optional
        Coordinate pair around which to center the layout. Default is the
        origin.
    dim : `int`
        Dimension of layout. If `dim` > 2, the remaining dimensions are set to
        zero in the returned positions.

    Returns
    -------
    pos : mapping
        A mapping of positions keyed by node.

    Examples
    --------
    >>> import pandas as pd
    >>> from pyrrole import ChemicalSystem
    >>> from pyrrole.drawing import diagram_layout
    >>> data = pd.DataFrame(
    ...     [{"name": "Separated_Reactants", "freeenergy": 0.},
    ...      {"name": "mlC1", "freeenergy": -5.4},
    ...      {"name": "mlC2", "freeenergy": -15.6},
    ...      {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
    ...      {"name": "mCARB1", "freeenergy": -9.7},
    ...      {"name": "mCARB2", "freeenergy": -19.8},
    ...      {"name": "mCARBX", "freeenergy": 20}]).set_index("name")
    >>> system = ChemicalSystem(
    ...     ["Separated_Reactants -> mlC1 -> mTS1",
    ...      "Separated_Reactants -> mlC2 -> mTS1",
    ...      "mCARB2 <- mTS1 -> mCARB1",
    ...      "Separated_Reactants -> mCARBX"], data)
    >>> digraph = system.to_digraph()
    >>> layout = diagram_layout(digraph)
    >>> layout['mCARB2']
    array([  3. , -19.8])

    Passing ``scale=1`` means scaling positions to ``(-1, 1)`` in all axes:

    >>> layout = diagram_layout(digraph, scale=1)
    >>> layout['mTS1'][1] <= 1.
    True

    """
    # TODO: private function of packages should not be used.
    graph, center = _nx.drawing.layout._process_params(graph, center, dim)

    num_nodes = len(graph)
    if num_nodes == 0:
        return {}
    elif num_nodes == 1:
        return {_nx.utils.arbitrary_element(graph): center}

    if sources is None:
        sources = graph.nodes()
    if targets is None:
        targets = graph.nodes()
    simple_paths = [path for source in set(sources) for target in set(targets)
                    for path in _nx.all_simple_paths(graph, source, target)]

    if pos is None:
        pos = tower_layout(graph, height=height, scale=None, center=center,
                           dim=dim)

    for path in simple_paths:
        for n, step in enumerate(path):
            if pos[step][0] < n:
                pos[step][0] = n

    if scale is not None:
        pos_arr = _np.array([pos[node] for node in graph])
        pos_arr = _nx.drawing.layout.rescale_layout(pos_arr,
                                                    scale=scale) + center
        pos = dict(zip(graph, pos_arr))

    # TODO: make test
    return pos


def draw_diagram_nodes(graph, pos=None, nodelist=None, node_size=.7,
                       node_color='k', style='solid', alpha=1.0, cmap=None,
                       vmin=None, vmax=None, ax=None, label=None):
    """
    Draw nodes of graph.

    This draws only the nodes of graph as horizontal lines at each
    ``y = pos[1]`` from ``x - node_size/2`` to ``x + node_size/2``, where
    ``x = pos[0]``.

    Parameters
    ----------
    graph : `networkx.Graph`
        A NetworkX graph.
    pos : mapping, optional
        A mapping with nodes as keys and positions as values. Positions should
        be sequences of length 2. If not specified (default), a diagram layout
        positioning will be computed. See `networkx.layout` and
        `pyrrole.drawing` for functions that compute node positions.
    nodelist : `list`, optional
        Draw only specified nodes (default is ``graph.nodes()``).
    node_size : scalar or array
        Size of nodes (default is ``.7``). If an array is specified it must be
        the same length as nodelist.
    node_color : color `str`, or array of `float`
        Node color. Can be a single color format `str` (default is ``'k'``), or
        a  sequence of colors with the same length as nodelist. If numeric
        values are specified they will be mapped to colors using the `cmap` and
        `vmin`, `vmax` parameters. See `matplotlib.hlines` for more details.
    style : `str` (``'solid'``, ``'dashed'``, ``'dotted'``, ``'dashdot'``)
        Edge line style (default is ``'solid'``). See `matplotlib.hlines` for
        more details.
    alpha : `float` or array of `float`, optional
        The node transparency. This can be a single alpha value (default is
        ``'1.0'``), in which case it will be applied to all the nodes of color.
        Otherwise, if it is an array, the elements of alpha will be applied to
        the colors in order (cycling through alpha multiple times if
        necessary).
    cmap : Matplotlib colormap, optional
        Colormap name or Colormap instance for mapping intensities of nodes.
    vmin : `float`, optional
        Minimum for node colormap scaling.
    vmax : `float`, optional
        Maximum for node colormap scaling.
    ax : `matplotlib.axes.Axes`, optional
        Draw the graph in the specified Matplotlib axes.
    label : `str`,  optional
        Label for legend.

    Returns
    -------
    `matplotlib.collections.LineCollection`
        `LineCollection` of the nodes.

    Examples
    --------
    >>> import pandas as pd
    >>> from pyrrole import ChemicalSystem
    >>> from pyrrole.drawing import draw_diagram_nodes
    >>> data = pd.DataFrame(
    ...     [{"name": "Separated_Reactants", "freeenergy": 0.},
    ...      {"name": "mlC1", "freeenergy": -5.4},
    ...      {"name": "mlC2", "freeenergy": -15.6},
    ...      {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
    ...      {"name": "mCARB1", "freeenergy": -9.7},
    ...      {"name": "mCARB2", "freeenergy": -19.8},
    ...      {"name": "mCARBX", "freeenergy": 20}]).set_index("name")
    >>> system = ChemicalSystem(
    ...     ["Separated_Reactants -> mlC1 -> mTS1",
    ...      "Separated_Reactants -> mlC2 -> mTS1",
    ...      "mCARB2 <- mTS1 -> mCARB1",
    ...      "Separated_Reactants -> mCARBX"], data)
    >>> digraph = system.to_digraph()
    >>> nodes = draw_diagram_nodes(digraph)

    """
    if ax is None:
        ax = _plt.gca()

    if nodelist is None:
        nodelist = list(graph.nodes())

    if not nodelist or len(nodelist) == 0:  # empty nodelist, no drawing
        return None

    if pos is None:
        pos = diagram_layout(graph)

    try:
        xy = _np.asarray([pos[v] for v in nodelist])
    except KeyError as e:
        raise _nx.NetworkXError('Node {} has no position.'.format(e))
    except ValueError:
        raise _nx.NetworkXError('Bad value in node positions.')

    if isinstance(alpha, _collections.Iterable):
        node_color = _nx.drawing.apply_alpha(node_color, alpha, nodelist, cmap,
                                             vmin, vmax)
        alpha = None

    node_collection = ax.hlines(xy[:, 1],
                                xy[:, 0] - node_size/2.,
                                xy[:, 0] + node_size/2.,
                                colors=node_color,
                                linestyles=style,
                                label=label,
                                cmap=cmap)

    node_collection.set_zorder(2)
    return node_collection


def draw_diagram_edges(graph, pos=None, edgelist=None, width=1.0,
                       edge_color='k', style='dashed', alpha=1.0,
                       edge_cmap=None, edge_vmin=None, edge_vmax=None, ax=None,
                       label=None, nodelist=None, node_size=.7):
    """
    Draw edges of graph.

    This draws only the edges of a graph.

    Parameters
    ----------
    graph : `networkx.Graph`
        A NetworkX graph.
    pos : mapping, optional
        A mapping with nodes as keys and positions as values. Positions should
        be sequences of length 2. If not specified (default), a diagram layout
        positioning will be computed. See `networkx.layout` and
        `pyrrole.drawing` for functions that compute node positions.
    edgelist : collection of edge `tuple`
        Draw only specified edges (default is ``graph.edges()``).
    width : `float`, or array of `float`
        Line width of edges (default is ``1.0``).
    edge_color : color `str`, or array of `float`
        Edge color. Can be a single color format `str` (default is ``'r'``),
        or a sequence of colors with the same length as edgelist. If numeric
        values are specified they will be mapped to colors using the
        `edge_cmap` and `edge_vmin`, `edge_vmax` parameters.
    style : `str` (``'solid'``, ``'dashed'``, ``'dotted'``, ``'dashdot'``)
        Edge line style (default is ``'dashed'``). See `matplotlib.hlines` for
        more details.
    alpha : `float`, optional
        The edge transparency (default is ``1.0``).
    edge_cmap : Matplotlib colormap, optional
        Colormap for mapping intensities of edges.
    edge_vmin : `float`, optional
        Minimum for edge colormap scaling.
    edge_vmax : `float`, optional
        Maximum for edge colormap scaling.
    ax : `matplotlib.axes.Axes`, optional
        Draw the graph in the specified Matplotlib axes.
    label : `str`,  optional
        Label for legend.
    nodelist : `list`, optional
        Draw only specified nodes (default is ``graph.nodes()``).
    node_size : scalar or array
        Size of nodes (default is ``.7``). If an array is specified it must be
        the same length as nodelist.

    Returns
    -------
    `matplotlib.collections.LineCollection`
        `LineCollection` of the edges.

    Examples
    --------
    >>> import pandas as pd
    >>> from pyrrole import ChemicalSystem
    >>> from pyrrole.drawing import draw_diagram_edges
    >>> data = pd.DataFrame(
    ...     [{"name": "Separated_Reactants", "freeenergy": 0.},
    ...      {"name": "mlC1", "freeenergy": -5.4},
    ...      {"name": "mlC2", "freeenergy": -15.6},
    ...      {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
    ...      {"name": "mCARB1", "freeenergy": -9.7},
    ...      {"name": "mCARB2", "freeenergy": -19.8},
    ...      {"name": "mCARBX", "freeenergy": 20}]).set_index("name")
    >>> system = ChemicalSystem(
    ...     ["Separated_Reactants -> mlC1 -> mTS1",
    ...      "Separated_Reactants -> mlC2 -> mTS1",
    ...      "mCARB2 <- mTS1 -> mCARB1",
    ...      "Separated_Reactants -> mCARBX"], data)
    >>> digraph = system.to_digraph()
    >>> edges = draw_diagram_edges(digraph)

    """
    if ax is None:
        ax = _plt.gca()

    if edgelist is None:
        edgelist = list(graph.edges())

    if not edgelist or len(edgelist) == 0:  # no edges!
        return None

    if nodelist is None:
        nodelist = list(graph.nodes())

    if pos is None:
        pos = diagram_layout(graph)

    try:
        # set edge positions
        edge_pos = _np.asarray([(pos[e[0]] + node_size/2.,
                                 pos[e[1]] - node_size/2.) for e in edgelist])
    except KeyError as e:
        raise _nx.NetworkXError('Node {} has no position.'.format(e))
    except ValueError:
        raise _nx.NetworkXError('Bad value in node positions.')

    if not _cb.iterable(width):
        lw = (width,)
    else:
        lw = width

    if not isinstance(edge_color, str) \
            and _cb.iterable(edge_color) \
            and len(edge_color) == len(edge_pos):
        if _np.alltrue([isinstance(c, str) for c in edge_color]):
            # (should check ALL elements)
            # list of color letters such as ['k','r','k',...]
            edge_colors = tuple([_colorConverter.to_rgba(c, alpha)
                                 for c in edge_color])
        elif _np.alltrue([not isinstance(c, str) for c in edge_color]):
            # If color specs are given as (rgb) or (rgba) tuples, we're OK
            if _np.alltrue([_cb.iterable(c) and len(c) in (3, 4)
                            for c in edge_color]):
                edge_colors = tuple(edge_color)
            else:
                # numbers (which are going to be mapped with a colormap)
                edge_colors = None
        else:
            raise ValueError('edge_color must contain color names or numbers')
    else:
        if isinstance(edge_color, str) or len(edge_color) == 1:
            edge_colors = (_colorConverter.to_rgba(edge_color, alpha), )
        else:
            raise ValueError('edge_color must be a color or list of one color '
                             ' per edge')

    edge_collection = _LineCollection(edge_pos,
                                      colors=edge_colors,
                                      linewidths=lw,
                                      antialiaseds=(1,),
                                      linestyle=style,
                                      transOffset=ax.transData)

    edge_collection.set_zorder(1)  # edges go behind nodes
    edge_collection.set_label(label)
    ax.add_collection(edge_collection)

    if _cb.is_numlike(alpha):
        edge_collection.set_alpha(alpha)

    if edge_colors is None:
        if edge_cmap is not None:
            assert(isinstance(edge_cmap, _Colormap))
        edge_collection.set_array(_np.asarray(edge_color))
        edge_collection.set_cmap(edge_cmap)
        if edge_vmin is not None or edge_vmax is not None:
            edge_collection.set_clim(edge_vmin, edge_vmax)
        else:
            edge_collection.autoscale()

    ax.autoscale_view()
    return edge_collection


def draw_diagram_labels(graph, pos=None, labels=None, font_size=12,
                        font_color='k', font_family='sans-serif',
                        font_weight='normal', alpha=1.0, bbox=None, ax=None,
                        offset=None, **kwds):
    """
    Draw node labels of graph.

    This draws only the node labels of a graph.

    Parameters
    ----------
    graph : `networkx.Graph`
        A NetworkX graph.
    pos : mapping, optional
        A mapping with nodes as keys and positions as values. Positions should
        be sequences of length 2. If not specified (default), a diagram layout
        positioning will be computed. See `networkx.layout` and
        `pyrrole.drawing` for functions that compute node positions.
    labels : mapping, optional
        Node labels in a mapping keyed by node of text labels.
    font_size : `int`, optional
       Font size for text labels (default is ``12``).
    font_color : `str`, optional
       Font color `str` (default is ``'k'``, i.e., black).
    font_family : `str`, optional
       Font family (default is ``'sans-serif'``).
    font_weight : `str`, optional
       Font weight (default is ``'normal'``).
    alpha : `float`, optional
        The text transparency (default is ``1.0``).
    ax : `matplotlib.axes.Axes`, optional
        Draw the graph in the specified Matplotlib axes.
    offset : array-like or `str`, optional
        Label positions are summed to this before drawing. Defaults to zero
        vector. If `str`, can be either ``'above'`` (equivalent to
        ``(0, 1.5)``) or ``'below'`` (equivalent to ``(0, -1.5)``).

    Returns
    -------
    mapping
        Mapping of labels keyed on the nodes.

    Examples
    --------
    >>> import pandas as pd
    >>> from pyrrole import ChemicalSystem
    >>> from pyrrole.drawing import draw_diagram_labels
    >>> data = pd.DataFrame(
    ...     [{"name": "Separated_Reactants", "freeenergy": 0.},
    ...      {"name": "mlC1", "freeenergy": -5.4},
    ...      {"name": "mlC2", "freeenergy": -15.6},
    ...      {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
    ...      {"name": "mCARB1", "freeenergy": -9.7},
    ...      {"name": "mCARB2", "freeenergy": -19.8},
    ...      {"name": "mCARBX", "freeenergy": 20}]).set_index("name")
    >>> system = ChemicalSystem(
    ...     ["Separated_Reactants -> mlC1 -> mTS1",
    ...      "Separated_Reactants -> mlC2 -> mTS1",
    ...      "mCARB2 <- mTS1 -> mCARB1",
    ...      "Separated_Reactants -> mCARBX"], data)
    >>> digraph = system.to_digraph()
    >>> edges = draw_diagram_labels(digraph, font_color='blue',
    ...                             offset="below")
    >>> labels = {k: "{:g}".format(v)
    ...           for k, v in digraph.nodes(data='freeenergy')}
    >>> edges = draw_diagram_labels(digraph, labels=labels,
    ...                             offset="above")

    """
    if ax is None:
        ax = _plt.gca()

    if labels is None:
        labels = dict((n, n) for n in graph.nodes())

    if pos is None:
        pos = diagram_layout(graph)

    if offset is None:
        offset = _np.array([0., 0.])
    elif offset == "above":
        offset = _np.array([0., 1.5])
    elif offset == "below":
        offset = _np.array([0., -1.5])

    # set optional alignment
    horizontalalignment = kwds.get('horizontalalignment', 'center')
    verticalalignment = kwds.get('verticalalignment', 'center')

    text_items = {}  # there is no text collection so we'll fake one
    for n, label in labels.items():
        (x, y) = _np.asanyarray(pos[n]) + _np.asanyarray(offset)
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same
        t = ax.text(x, y, label,
                    size=font_size,
                    color=font_color,
                    family=font_family,
                    weight=font_weight,
                    alpha=alpha,
                    horizontalalignment=horizontalalignment,
                    verticalalignment=verticalalignment,
                    transform=ax.transData,
                    bbox=bbox,
                    clip_on=True)
        text_items[n] = t

    return text_items


def draw_diagram(graph, pos=None, with_labels=True, offset=None, **kwds):
    """
    Draw a diagram for graph using Matplotlib.

    Draw graph as a simple energy diagram with Matplotlib with options for node
    positions, labeling, titles, and many other drawing features. See examples
    below.

    Parameters
    ----------
    graph : `networkx.Graph`
        A NetworkX graph.
    pos : mapping, optional
        A mapping with nodes as keys and positions as values. Positions should
        be sequences of length 2. If not specified (default) a diagram layout
        positioning will be computed. See `networkx.drawing.layout` and
        `pyrrole.drawing` for functions that compute node positions.
    with_labels : `bool`, optional
       Set to `True` (default) to draw labels on the nodes.
    offset : array-like or `str`, optional
        Label positions are summed to this before drawing. Defaults to
        ``'below'``. See `draw_diagram_labels` for more.

    Notes
    -----
    Further keywords are passed to `draw_diagram_nodes` and
    `draw_diagram_edges`. If `pos` is `None`, `diagram_layout` is also called
    and have keywords passed as well. The same happens with
    `draw_diagram_labels` if `with_labels` is `True`.

    Examples
    --------
    >>> import pandas as pd
    >>> from pyrrole import ChemicalSystem
    >>> from pyrrole.drawing import draw_diagram, draw_diagram_labels
    >>> data = pd.DataFrame(
    ...     [{"name": "Separated_Reactants", "freeenergy": 0.},
    ...      {"name": "mlC1", "freeenergy": -5.4},
    ...      {"name": "mlC2", "freeenergy": -15.6},
    ...      {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
    ...      {"name": "mCARB1", "freeenergy": -9.7},
    ...      {"name": "mCARB2", "freeenergy": -19.8},
    ...      {"name": "mCARBX", "freeenergy": 20}]).set_index("name")
    >>> system = ChemicalSystem(
    ...     ["Separated_Reactants -> mlC1 -> mTS1",
    ...      "Separated_Reactants -> mlC2 -> mTS1",
    ...      "mCARB2 <- mTS1 -> mCARB1",
    ...      "Separated_Reactants -> mCARBX"], data)
    >>> digraph = system.to_digraph()
    >>> draw_diagram(digraph)
    >>> labels = {k: "{:g}".format(v)
    ...           for k, v in digraph.nodes(data='freeenergy')}
    >>> edges = draw_diagram_labels(digraph, labels=labels,
    ...                             offset="above")

    """
    if pos is None:
        pos = diagram_layout(graph, **kwds)  # default to diagram layout

    node_collection = draw_diagram_nodes(graph, pos, **kwds)  # noqa
    edge_collection = draw_diagram_edges(graph, pos, **kwds)  # noqa
    if with_labels:
        if offset is None:
            # TODO: This changes the default behaviour of draw_diagram_labels.
            offset = "below"
        draw_diagram_labels(graph, pos, offset=offset, **kwds)
    _plt.draw_if_interactive()
