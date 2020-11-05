#!/usr/bin/python
import networkx as nx
import graph_tool as gt

import graph_tool.draw
import graph_tool.topology

import math

import time as tm # to see what takes a long time.


DEBUG = False;

OPEN = False;
CLOSED = True;

class Interval:

    def __init__(self, left_boundary, left, right, right_boundary):
        self.left_boundary = left_boundary;
        self.left = left;
        self.right = right;
        self.right_boundary = right_boundary;

    # returns interval we can use in the split graph generated in the
    # reeb graph class.
    # TODO: must modify if allowing non integer time intervals in rg.
    def get_split_interval(self):
        split = [];
        if self.left_boundary == OPEN:
            split.append(self.left + .5);
        else:
            split.append(self.left);


        if self.right_boundary == OPEN:
            split.append(self.right - .5);
        else:
            split.append(self.right);

        return Interval(CLOSED, split[0], split[1], CLOSED);

    def minus(self):
        left_boundary = self.left_boundary;
        left_value = self.left;

        if self.left_boundary == OPEN:
            left_boundary = CLOSED;
        else:
            left_value = left_value - 1;
            left_boundary = CLOSED;

        return Interval(left_boundary, left_value, self.right, self.right_boundary);

    def plus(self):
        right_boundary = self.right_boundary;
        right_value = self.right;

        if self.right_boundary == OPEN:
            right_boundary = CLOSED;
        else:
            right_value = right_value + 1;
            right_boundary = OPEN;

        return Interval(self.left_boundary, self.left, right_value, right_boundary);

    def plus_minus(self):
        left_boundary = self.left_boundary;
        left_value = self.left;

        if self.left_boundary == OPEN:
            left_boundary = CLOSED;
        else:
            left_value = left_value - 1;
            left_boundary = CLOSED;

        right_boundary = self.right_boundary;
        right_value = self.right;

        if self.right_boundary == OPEN:
            right_boundary = CLOSED;
        else:
            right_value = right_value + 1;
            right_boundary = OPEN;


        return Interval(left_boundary, left_value, right_value, right_boundary);

    def __str__(self):
        representation = "";
        if self.left_boundary == OPEN:
            representation += "(";
        else:
            representation += "[";

        representation += str(self.left) + "," + str(self.right)

        if self.right_boundary == OPEN:
            representation += ")";
        else:
            representation += "]";

        return representation;

class Node:

    def __init__(self, time, label):
        self.time = time;
        self.label = label;

    def get_time(self):
        return self.time;

    def get_label(self):
        return self.label;

    def copy(self):
        return Node(self.time, self.label);

    # assumes this is an interval object
    def in_closed_interval(self, interval):
        # assuming closed so I don't have to add too many checks;
        # need this to use specifically with closed intervals for now.
        is_within = False;

        if interval.left <= self.time and self.time <= interval.right:
            is_within = True;

        return is_within;

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.time == other.time and self.label == other.label;
        else:
            return False;

    def __str__(self):
        return "Node({}, {})".format(self.time, self.label);

class Edge:

    def __init__(self, source, target):
        self.source = source; # source vertex
        self.target = target; # target vertex

    def get_source(self):
        return self.source;

    def get_target(self):
        return self.target;

    def copy(self):
        return Edge(self.source.copy(), self.target.copy());

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.source == other.source and self.target == other.target;
        else:
            return False;

    def __str__(self):
        return "Edge({}, {})".format(self.source, self.target);

class Graph:

    def __init__(self, nodes, edges):
        self.nodes = nodes;
        self.edges = edges;

        # throw error if nodes not integers
        # throw error if edges not directed

    def __str__(self):
        str = ""
        for node in self.nodes:
            str += "{}\n".format(node);

        for edge in self.edges:
            str += "{}\n".format(edge);

        return str;

    def get_nodes(self):
        return self.nodes;

    def get_edges(self):
        return self.edges;

    def add_node(self, node):
        return None;

    def add_edge(self, edge):
        return None;

    def contains_edge(self, edge_source_time, edge_target_time):
        # edge_exists = False;
        for edge in self.edges:
            if edge.get_source().get_time() == edge_source_time and edge.get_target().get_time() == edge_target_time:
                return True;
        return False;


    def get_graph_tool(self):
        # TODO : need to save node data internally to graph.
        # this is needed to reed back the component labeling
        g = gt.Graph();

        time = g.new_vertex_property("object");
        g.vertex_properties["time"] = time;

        label = g.new_vertex_property("object");
        g.vertex_properties["label"] = label;

        # g.list_properties()

        # TODO: might add vertex start end data too, but not necessary atm.
        # can be recovered from the vertex data.

        # TODO: consider deleting the dictionaries as they're not being used.

        nodes = {};
        for node in self.nodes:
            vertex = g.add_vertex()
            nodes[(node.get_time(), node.get_label())] = vertex;

            # store attributes internally to graph
            g.vp.time[vertex] = node.get_time();
            g.vp.label[vertex] = node.get_label();

        edges = {};
        for edge in self.edges:
            source = edge.get_source()
            source_time = source.get_time();
            source_label = source.get_label();

            target = edge.get_target();
            target_time = target.get_time();
            target_label = target.get_label();

            edges[(source_time, source_label, target_time, target_label)] = g.add_edge(nodes[(source_time, source_label)], nodes[(target_time, target_label)]);

        return g;

    # assumes title doesn't have .pdf at end
    def generate_pdf(self, title):
        filename = str(title) + ".pdf";
        g = self.get_graph_tool();
        graph_tool.draw.graph_draw(g, vertex_text=g.vertex_index, output=filename);

# assume critical points at integers
class ReebGraph:

    # creates ReebGraph with empty graph
    # assumes integer labels
    def __init__(self, nodes, edges):
        # self.directedGraph = nx.Graph();
        # self.nodes = nodes;
        # self.edges = edges;

        # save critical points for barcode calculation :
        critical_times = []
        for node in nodes:
            if node.get_time() not in critical_times:
                critical_times.append(node.get_time());

        self.critical_times = critical_times;

        self.reeb_graph = Graph(nodes, edges);
        if DEBUG: print(self.reeb_graph);

        self.partition = ReebGraph.partition(self.reeb_graph);
        if DEBUG: print(self.partition);

        self.split = ReebGraph.split(self.partition);
        if DEBUG: print(self.split)

    # def __str__(self):

    def get_next_node_label(nodes, t):
        label = -1;
        for node in nodes:
            node_time = node.get_time();
            node_label = node.get_label();

            if node_time == t and label < node_label:
                label = node_label;

        return label + 1;

    # takes Reeb Graph and breaks edges to be dimention one
    def partition(graph):

        nodes = graph.nodes;
        edges = graph.edges;

        nodes_p = [];
        edges_p = [];

        for node in nodes:
            # if DEBUG: print(node);
            nodes_p.append(node);

        for edge in edges:
            edge_source = edge.get_source();
            edge_target = edge.get_target();

            start_time = edge_source.get_time();
            end_time = edge_target.get_time();

            if (start_time + 1) == end_time:
                # if DEBUG: print("{} and {} are one unit apart".format(edge_source, edge_target));
                edges_p.append(edge);
            else:
                # delete edge and add new ones based on start and end time:
                # if DEBUG: print("{} bigger than one unit. Must split.".format(edge));

                last_node = edge_source;
                for t in range(start_time + 1, end_time):
                    # if DEBUG: print("Add node at time {}. Add edge before".format(t));

                    new_node = Node(t, ReebGraph.get_next_node_label(nodes_p, t))

                    # add new node
                    nodes_p.append(new_node);

                    # add new edges to this node
                    edges_p.append(Edge(last_node, new_node));
                    last_node = new_node;

                # add last edge
                edges_p.append(Edge(last_node, edge_target));

            # print(edge)

        # for node in nodes_p:
        #     print(node)
        # for edge in edges_p:
        #     print(edge)
        return Graph(nodes_p, edges_p);

    def split(graph):

        # print("here")
        # print(graph);

        nodes = graph.nodes;
        edges = graph.edges;

        nodes_s = []
        edges_s = []

        for node in nodes:
            # if DEBUG: print(node);
            nodes_s.append(node);

        for edge in edges:
            edge_source = edge.get_source();
            edge_target = edge.get_target();

            start_time = edge_source.get_time();
            end_time = edge_target.get_time();

            # add intermediate node
            intermediate_time = start_time + (end_time - start_time) / 2
            intermediate_node = Node(intermediate_time, ReebGraph.get_next_node_label(nodes_s, intermediate_time)); # might have labeling problems this way?
            nodes_s.append(intermediate_node);

            edges_s.append(Edge(edge_source, intermediate_node));
            edges_s.append(Edge(intermediate_node, edge_target));


        return Graph(nodes_s, edges_s);


    def get_in_degree(self, node):
        in_degree = 0;

        for edge in self.reeb_graph.edges:
            if edge.get_target() == node:
                in_degree += 1

        return in_degree;

    def get_out_degree(self, node):
        out_degree = 0;

        for edge in self.reeb_graph.edges:
            if edge.get_source() == node:
                out_degree += 1

        return out_degree;

    # TODO : consider deleting; not in use.
    def get_vertices_at_time(self, time):
        vertex_count = 0;

        # is time a critical point?
        if math.floor(time) == time:
            if DEBUG: print("{} is integer".format(time));
        else:
            time = math.floor(time) + 0.5;

        for node in self.split.nodes:
            if node.get_time() == time:
                vertex_count += 1;


        return vertex_count;

    def get_edges_between_times(self, start_time, end_time):
        edges = [];

        return edges;

    def get_subgraph(self, interval):
        split_interval = interval.get_split_interval();
        split_start = split_interval.left;
        split_end = split_interval.right;

        subgraph_nodes = [];
        subgraph_edges = [];

        # code asumes critical points are the integers, will need to keep track
        # of them in the future if want arbitary ones.

        # generate all nodes in between the interval
        for node in self.split.nodes:

            is_integer = node.time == math.floor(node.time);
            if node.in_closed_interval(split_interval):
                # if DEBUG: print(node);
                subgraph_nodes.append(node);

        # generate all edges in between the interval
        for edge in self.split.edges:

            if edge.source in subgraph_nodes and edge.target in subgraph_nodes:
                # if DEBUG: print(edge);
                subgraph_edges.append(edge);

        return Graph(subgraph_nodes, subgraph_edges);


    # assumes the time data of the vertices in a connected grapha are given.
    # edges unnecessary.
    # assumes interval is closed.
    def is_connected_subgraph_full(self, subgraph_node_times, closed_interval):

        exists_node_at_interval_start = False;
        exists_node_at_interval_end = False;

        for time in subgraph_node_times:
            if time == closed_interval.left:
                exists_node_at_interval_start = True;
            if time == closed_interval.right:
                exists_node_at_interval_end = True;

        return exists_node_at_interval_start and exists_node_at_interval_end;

    # interval assumed to be interval object
    def full(self, interval):
        # preprocesses inputs
        subgraph = self.get_subgraph(interval);

        split_interval = interval.get_split_interval();

        # subgraph.generate_pdf("test");
        # print(subgraph)
        # print(interval);
        # print(split_interval);

        # we pass this subgraph into the graph-tool package to calculate
        # connected components for us.
        g = subgraph.get_graph_tool();
        components = graph_tool.topology.label_components(g, directed=False)[0].a

        if DEBUG: print(components);

        # code below assumes component labeling starts from 0 and increases.
        connected_subgraphs = []
        for i in range(0, len(set(components))):
            connected_subgraphs.append([]);

        for i, component in enumerate(components):

            node_time = g.vp.time[g.vertex(i)];
            node_label = g.vp.label[g.vertex(i)];

            connected_subgraphs[component].append(node_time);

            if DEBUG: print("Node:({}, {}) is in component {}".format(node_time, node_label, g.vertex(component)));

        # donw with graph-tool object

        # print(connected_subgraphs)

        full_counter = 0;
        for subgraph in connected_subgraphs:
            # print(subgraph);
            if self.is_connected_subgraph_full(subgraph, split_interval):
                full_counter += 1;

        # print(full_counter);
        return full_counter;

    def multiplicity(self, interval):
        # print(interval);
        # print(interval.plus())
        # print(interval.minus())
        # print(interval.plus_minus())

        multiplicity = self.full(interval) + self.full(interval.plus_minus()) \
            - self.full(interval.plus()) \
            - self.full(interval.minus());

        # if multiplicity != 0:
        #     i = interval;
        #     print("{} has {} full components.".format(i, self.full(i)));
        #     i = interval.plus();
        #     print("{} has {} full components.".format(i, self.full(i)));
        #     i = interval.minus();
        #     print("{} has {} full components.".format(i, self.full(i)));
        #     i = interval.plus_minus();
        #     print("{} has {} full components.".format(i, self.full(i)));

        return multiplicity;

    def barcode(self):
        # this should return a dictionary with all intervals in the lifetime
        # of the reeb graph. The points in the intervals of interest are
        # all possible critical points, so, for all nodes in the partitioned
        # graph.

        intervals = [];

        for i in range(0, len(self.critical_times)):
            for j in range(i + 1, len(self.critical_times)):
                # print((i, j))

                left = self.critical_times[i];
                right =  self.critical_times[j];

                # check if this corresponds to an edge.
                # if self.reeb_graph.contains_edge(left, right):
                intervals.append(Interval(CLOSED, left, right, CLOSED));
                intervals.append(Interval(OPEN, left, right, CLOSED));
                intervals.append(Interval(CLOSED, left, right, OPEN));
                intervals.append(Interval(OPEN, left, right, OPEN));

        barcode = {};
        for interval in intervals:
            full = self.multiplicity(interval);
            if full != 0:
                print("{} has multiplicity {}".format(interval, full));

        return None;


nodes = [
    Node(0, 0),
    Node(1, 0),
    Node(2, 0),
    Node(3, 0),
    Node(4, 0),
    Node(5, 0)
    ];
edges = [
    Edge(Node(0, 0), Node(2, 0)),
    Edge(Node(1, 0), Node(2, 0)),
    Edge(Node(2, 0), Node(3, 0)),
    Edge(Node(3, 0), Node(4, 0)),
    Edge(Node(3, 0), Node(4, 0)),
    Edge(Node(4, 0), Node(5, 0))
]

time_start = tm.time();

rg = ReebGraph(nodes, edges);

# for time in range(0, 6):
#     print("time:{} | in:{} | out:{}".format(time, rg.get_in_degree(Node(time, 0)), rg.get_out_degree(Node(time, 0))));
#
# g = rg.reeb_graph.get_graph_tool();
# print(graph_tool.topology.label_components(g, directed=False)[0].a);
# graph_tool.draw.graph_draw(g, vertex_text=g.vertex_index, output="ReebGraph.pdf");
#
# g = rg.partition.get_graph_tool();
# graph_tool.draw.graph_draw(g, vertex_text=g.vertex_index, output="Partition.pdf");
#
# g = rg.split.get_graph_tool();
# graph_tool.draw.graph_draw(g, vertex_text=g.vertex_index, output="Split.pdf");

# test code for interval class
# interval = Interval(OPEN, 1, 2, CLOSED);
# print(interval);
# print(interval.plus())
# print(interval.minus())
# print(interval.plus_minus())
# print(interval.get_split_interval());

# print(rg.get_vertices_at_time(1));

subgraph = rg.get_subgraph(Interval(OPEN, 1, 4, OPEN));
# subgraph.generate_pdf(tm.time());
# print(subgraph);

time_end = tm.time();
print("runtime : " + str(time_end - time_start) + " seconds")

# rg to test full function calculations


nodes = [
    Node(1, 0),
    Node(2, 0),
    Node(3, 0),
    Node(4, 0)
    ];
edges = [
    Edge(Node(1, 0), Node(3, 0)),
    Edge(Node(2, 0), Node(3, 0)),
    Edge(Node(2, 0), Node(4, 0))
]
rg = ReebGraph(nodes, edges);
rg.barcode()
print()

time_start = tm.time();

# Bei Wang example
# nodes = [
#     Node(0, 0),
#     Node(1, 0),
#     Node(2, 0),
#     Node(3, 0),
#     Node(4, 0),
#     Node(5, 0),
#     Node(6, 0),
#     Node(7, 0)
#     ];
# edges = [
#     Edge(Node(0, 0), Node(2, 0)),
#     Edge(Node(1, 0), Node(2, 0)),
#     Edge(Node(2, 0), Node(3, 0)),
#     Edge(Node(3, 0), Node(5, 0)),
#     Edge(Node(3, 0), Node(4, 0)),
#     Edge(Node(4, 0), Node(6, 0)),
#     Edge(Node(4, 0), Node(7, 0))
# ];


# Justin's example
nodes = [
    Node(0, 0),
    Node(1, 0),
    Node(2, 0),
    Node(3, 0),
    Node(4, 0),
    Node(5, 0)
    ];
edges = [
    Edge(Node(0, 0), Node(3, 0)),
    Edge(Node(1, 0), Node(4, 0)),
    Edge(Node(2, 0), Node(3, 0)),
    Edge(Node(3, 0), Node(4, 0)),
    Edge(Node(4, 0), Node(5, 0))
];
rg = ReebGraph(nodes, edges);
rg.barcode()
print()


nodes = [
    Node(0, 0),
    Node(1, 0),
    Node(2, 0),
    Node(3, 0),
    Node(4, 0),
    Node(5, 0)
    ];
edges = [
    Edge(Node(0, 0), Node(4, 0)),
    Edge(Node(1, 0), Node(3, 0)),
    Edge(Node(2, 0), Node(3, 0)),
    Edge(Node(3, 0), Node(4, 0)),
    Edge(Node(4, 0), Node(5, 0))
];
rg = ReebGraph(nodes, edges);
# rg.reeb_graph.generate_pdf("ReebGraph" + str(tm.time()));
# rg.partition.generate_pdf("Partition" + str(tm.time()));
# rg.split.generate_pdf("Split" + str(tm.time()));

# full = rg.full(Interval(CLOSED, 2, 3, OPEN));
# print(full)
# print(full);

# subgraph = rg.get_subgraph(Interval(CLOSED, 2, 3, CLOSED));
# print(subgraph);
# rg.reeb_graph.generate_pdf("ReebGraph" + str(tm.time()));
# rg.partition.generate_pdf("Partition" + str(tm.time()));
# rg.split.generate_pdf("Split" + str(tm.time()));
# print(rg.full(Interval(OPEN, 2, 4, OPEN)))
rg.barcode()

time_end = tm.time();
print("runtime : " + str(time_end - time_start) + " seconds")
