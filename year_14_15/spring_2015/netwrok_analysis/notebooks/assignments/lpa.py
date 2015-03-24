#! /usr/bin/python
"""Finds the communities via the label propagation algorithm."""

import networkx as nx
import os
import re
import sys
from optparse import OptionParser

__all__ = [ "semisynchronous_prec_max" , "print_communities" ]

USAGE = "%prog [OPTION ...] NETWORK_DESCRIPTION"
USE_DESCRIPTION = "The NETWORK_DESCRIPTION file must be an adjacency list."
parser = OptionParser( USAGE , description=__doc__ + " " + USE_DESCRIPTION )
parser.set_defaults( verbosity=None )
parser.add_option( "-a" , "--adj",
    action="store_const",
    const="adj",
    dest="reader_format",
    help="Indicates that the description file is in adjacency list format."
)
parser.add_option( "-g" , "--gml",
    action="store_const",
    const="gml",
    dest="reader_format",
    help="Indicates that the description file is in GML format."
)
parser.add_option( "-l" , "--log-verbosity",
    action="append",
    metavar="FILE",
    dest="verbosity",
    help="Verbosely inform the user of the ongoings of the computation "
         "logging to FILE. Last indicated verbosity level trumps all."
)
parser.add_option( "-v" , "--verbose",
    action="append_const",
    const="stdout",
    dest="verbosity",
    help="Verbosely inform the user of the ongoings of the computation to "
         "stdout. Last indicated verbosity level trumps all."
)

# mapping of file formats to networkx reading functions
network_reader = dict()
network_reader['adj'] = nx.read_adjlist
network_reader['gml'] = nx.read_gml

verbose_stream = None

def main():
    ( opts , args ) = parser.parse_args()
    __validate_opts_and_args( opts , args )
    __setup_verbosity( opts.verbosity )
    __inform_user( "parsing network ..." )
    network = network_reader[opts.reader_format]( args[0] )
    communities = semisynchronous_prec_max( network )
    print_communities( communities )

def semisynchronous_prec_max( network ):
    """Main functionality of the LPA algorithm.

       Factored out of the main method for external usage.
    """
    __inform_user( "coloring network ..." )
    coloring = __color_network( network )
    __inform_user( "uniquely labeling network ..." )
    labeling = __uniquely_label( network )
    nrounds = 1
    while not __labeling_complete( labeling , network ):
        nrounds += 1
        __update_labels( labeling , coloring , network )
    communities = __form_communities( labeling , network )
    return communities

def print_communities( communities ):
    """Prints the supplied dict() which maps labels to sets of nodes.

       Supplied dict() should be the results of an LPA community detection.
    """
    for l, nodes in communities.items():
        node_list = [ str(n) for n in nodes ]
        print re.sub( "[\[,\]]" , "" , "{0}: {1}".format(l,node_list) )

def __break_color_tie( current , labels ):
    """Uses Prec-Max tie-breaking to break the ties, as laid out in:
       'Community Detection via Semi-Synchronous Label Propagation Algorithms'
       Cordasco and Gargano, 2011

       If the labels set specified is empty than the current label is returned.

       Specified set of labels is assummed to be one of integers.
    """
    if len(labels) == 0 or current in labels:
        new_label = current
    else:
        # freakin' no built-in max function on sets!
        mx = -(sys.maxint) - 1
        for c in labels:
            mx = c if c > mx else mx
        new_label = mx
    return new_label

def __calculate_label_frequencies( node , labeling , network ):
    """Counts up the labels of the neighbors of the specified node. Returns a
       dictionary from the label to the frequency.
    """
    counts = dict()
    for q in network.neighbors( node ):
        qlabel = labeling[q]
        counts[qlabel] = 1 if qlabel not in counts else counts[qlabel] + 1
    return counts

def __color_network( net ):
    """Colors the network so that neighboring nodes all have distinct colors.
       Returns a dict of set of nodes and also a lookup of nodes to colors, in a
       tuple in that order.
    """
    coloring = dict() # color => set(node)
    lookup = dict() # node => color
    finally_colored_nodes = set()
    for n in net.nodes():
        color = 0 if n not in lookup else lookup[n]
        lookup[n] = color
        finally_colored_nodes.add( n )
        if color not in coloring: coloring[color] = set()
        coloring[color].add( n )
        for q in net.neighbors( n ):
            if q not in finally_colored_nodes: lookup[q] = lookup[n] + 1
    return coloring

def __form_communities( labeling , network ):
    """Determines the communities from the labels of the network, returning a
       dict of sets of nodes.
    """
    communities = dict() # label => set(nodes)
    for n in network.nodes():
        label = labeling[n]
        if label not in communities: communities[label] = set()
        communities[label].add( n )
    return communities

def __inform_user( msg ):
    """Writes supplied statement if verbosity is turned on."""
    if verbose_stream is not None: verbose_stream.write( msg + "\n" )

def __labeling_complete( labeling , network ):
    """Determines whether or not LPA is done. It is complete when all nodes have
       a label that is in the set of highest frequency labels amongst its
       neighbors.

       Nodes with no neighbors are themselves a community and are therefore
       labeled, hence the immediate if statement in the for loop.
    """
    for node in network:
        if len(network.neighbors(node)) != 0:
            counts = __calculate_label_frequencies( node , labeling , network )
            high_labels = __select_labels_of_highest_frequency( counts )
            if labeling[node] not in high_labels: return False
    return True

def __select_labels_of_highest_frequency( freqs ):
    """Finds all labels of maximum frequency. Specified freqs must be a mapping
       from label to frequency of that label.

       Returns a set.
    """
    labels = set()
    mx = -(sys.maxint) - 1
    for label, freq in freqs.items():
        if mx <= freq:
            if mx < freq:
                mx = freq
                labels.clear()
        labels.add( label )
    return labels

def __setup_verbosity( verbosity ):
    """Translates the command-line specified verbosity list (see program
       options) to the appropriate settings for verbosely explain the program's
       operations.

       'verbosity' is the list created by the commandline arguments.
    """
    global verbose_stream
    if verbosity is None:
        verbose_stream = None
    else:
        last_stream = verbosity[-1]
        verbose_stream = ( sys.stdout if last_stream=="stdout"
                                      else open(last_stream,"w") )

def __uniquely_label( network ):
    """Gives a unique label (integer) to each node in the network."""
    labeling = dict()
    for n, label in zip( network.nodes() , range(len(network))):
        labeling[n] = label
    return labeling

def __update_labels( labeling , coloring , network ):
    """Updates labels of every single node in the network."""
    for color, nodes in coloring.items():
        for n in nodes:
            __update_label( n , labeling , network )

def __update_label( node , labeling , network ):
    """Updates the label of a SINGLE node in the network."""
    counts = __calculate_label_frequencies( node , labeling , network )
    high_labels = __select_labels_of_highest_frequency( counts )
    if len(high_labels) == 1:
        labeling[node] = high_labels.pop()
    elif len(high_labels) > 1:
        labeling[node] = __break_color_tie( labeling[node] , high_labels )

def __validate_opts_and_args( opts , args ):
    """Makes sure that at least one file has been passed and that it exists.
       Also make sure that the version of the description file has been
       indicated via an option.
    """

    if opts.reader_format is None or len(args)<1 or len(args)>1:
        parser.print_help()
        sys.exit( 1 )
    if not os.path.exists( args[0] ): 
        sys.stderr.write( "File, {0}, does not exist.\n".format(args[0]) )
        parser.print_help()
        sys.exit( 1 )

if __name__ == "__main__":
    sys.exit( main() )
