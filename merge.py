
from simplicial import *
from graph_tool.all import *

DEBUG = False;

class MergeTree:

    def __init__(self):
        # return None
        if DEBUG : print("MergeTree : init")

        self.critical_times = []
        self.strands = []

class MergePresentation:

    def __init__(self, generators, relations):

        self.critical_times = set()
        self.generators = {} # generator_id : time
        self.relations = {} # relation_id : time

        # add generators as zero-simplices
        for i in range(len(generators)):
            # print(generators[i])
            self.critical_times.add(generators[i])
            # self.generators.append(str(i))
            self.generators[str(i)] = generators[i]

        # add relations as one-simplices
        for key, value in relations.items():
            u = str(key[0])
            v = str(key[1])

            # check valid filtration
            assert(self.generators[u] < value and self.generators[v] < value)

            self.relations[(u,v)] = value
            id_str = "{}:{}".format(u, v)

            self.critical_times.add(value)


    def __str__(self):
        pres_str = "M : {};\n".format(self.get_simplicial_complex().simplices())
        critical_str = "C_M := {};\n".format(str(self.critical_times));
        connected = "|pi_0(M)| = {}".format(self.count_connected_components())
        return pres_str + critical_str + connected
        # self.generators = []
        # self.relations = []

    def get_simplicial_complex(self):
        c = SimplicialComplex()

        # add generators as zero-simplices
        for i in range(len(generators)):
            # self.critical_times.add(generators[i])
            # self.generators.append(str(i))
            # self.generators[str(i)] = generators[i]
            c.addSimplex(id = str(i), attr = dict(time = generators[i]))

        # add relations as one-simplices
        for key, value in relations.items():
            u = str(key[0])
            v = str(key[1])

            # check valid filtration
            # assert(c[u]["time"] < value and c[v]["time"] < value)

            # self.relations.append((u, v))
            # self.relations[(u,v)] = value
            id_str = "{}:{}".format(u, v)
            c.addSimplex(fs = [u, v], id = id_str, attr = dict(time = value))

            # self.critical_times.add(value)
        return c

    def get_graph(self):
        g = Graph(directed=False)

        v_time = g.new_vertex_property("object");
        g.vertex_properties["time"] = v_time;

        e_time = g.new_edge_property("object")
        g.edge_properties["time"] = e_time

        # print(g.list_properties())

        vertices = {}
        for generator, time in self.generators.items():
            # print(generator)
            vertex = g.add_vertex()
            vertices[generator] = vertex
            g.vp.time[vertex] = time

        edges = {}
        for relation, time in self.relations.items():
            # print(relation)

            u = relation[0]
            v = relation[1]
            edge = g.add_edge(vertices[u], vertices[v])
            edges[relation] = edge
            g.ep.time[edge] = time

        return g

    def get_filtration(self):
        g = self.get_graph();

        for i in self.critical_times:
            vfilt=lambda v : g.vp.time[v] <= i
            efilt=lambda e: g.ep.time[e] <= i
            u = GraphView(g, vfilt=vfilt, efilt=efilt)
            # graph_draw(u, output="M[{}].pdf".format(i))


        return None

    def count_connected_components(self):
        g = self.get_graph()
        comp, hist = label_components(g)
        return len(set(comp.a))

    def is_connected(self):
        """
        Returns True if the merge forest is a tree.
        """
        return self.count_connected_components() == 1

generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
relations = {
    (0, 1) : 4,
    (1, 2) : 5
}

presentation = MergePresentation(generators, relations)
print(presentation)
# print(presentation.count_connected_components())
# print(presentation.is_connected())
g = presentation.get_graph()
presentation.get_filtration()

# filename = "test.pdf";
# graph_tool.draw.graph_draw(g, vertex_text=g.vp.label, output=filename);
# pos = sfdp_layout(g)
# deg = g.degree_property_map("in")
# graph_draw(g, pos=pos, vertex_size=5, output=filename)

# c = SimplicialComplex()
# # add a simplex with a generated name
# s1 = c.addSimplex()
#
# # add simplices whose names we want to specify
# s2 = c.addSimplex(id = 2)
# s3 = c.addSimplex(id = 3)
#
# l23 = c.addSimplex(fs = [ 2, 3 ])
#
# ss = c.simplices()
# print(ss)
