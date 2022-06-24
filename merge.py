
from poset_utils import *

from itertools import combinations
from simplicial import *
from graph_tool.all import *

DEBUG = False;

class Merge:

    def __init__(self, m, strands, merge_points):
        # return None
        # if DEBUG : print("Merge : init")
        self.dimension = m
        Rm = lambda c : R(m, c)

        self.set = {}

        # generate critical times set
        self.critical_times = set()
        for index, time in strands.items():
            self.critical_times.add(Rm(time))
        for indices, times in merge_points.items():
            for time in times:
                self.critical_times.add(Rm(time))

        # initialize internal strands
        self.strands = {}
        for index, time in strands.items():
            self.strands[index] = Rm(time)

        # initialize sets
        for time in self.critical_times:
            print(time)
            self.set[time] = downset_dict(self.strands, time)
            for element in self.set[time]:
                print("--{}".format(element))

        # # initialize sets
        # for index, time in strands.items():
        #     if Rm(time) not in self.set:
        #         self.set[Rm(time)] = []
        #     self.set[Rm(time)].append(index)
        #
        #     # self.strands[index] = Rm(time)
        # # self.strands = []
        # for key, value in self.set.items():
        #     print("{} - {}".format(key, value))
        #
        # # self.merge_points = {}
        # for indices, times in merge_points.items():
        #     # for index in indices:
        #     #     print(index)
        #     r_times = []
        #     for time in times:
        #         # print(time)
        #         r_times.append(Rm(time))
        #         self.critical_times.add(Rm(time))
        #     # self.merge_points[indices] = r_times

    def get_set(self, r):
        return None

    def get_structure_map(self, s, t):
        assert(s <= t)
        return None

    def get_structure_poset(self):
        poset = []
        for pair in combinations(self.critical_times, 2):
            x = pair[0]
            y = pair[1]
            if x <= y:
                poset.append([x, y])
            else:
                poset.append([y, x])
            # print(pair)
        return poset

    def __str__(self):
        text = ""
        components = []
        for time in self.critical_times:
            # print(time)
            components.append(str(time))
            # text += str(self.coordinate[i])
            # if i < self.dimension - 1:
            #     text += ","

        return "({})".format(", ".join(components))

m = 2
strands = {
    0 : [0, 0],
    1 : [0, 0],
    2 : [2, 2]
}

# generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
merge_points = {
    (0, 1, 2) : [[0, 1], [1, 0]]
}

merge = Merge(m, strands, merge_points)
poset = merge.get_structure_poset()
for pair in poset:
    print("{} <= {}".format(pair[0], pair[1]))
    # print(pair[0])
    # print(pair[1])
# print(merge)

exit()

class MergePresentation:

    def __init__(self, m, generators, relations):
        Rm = lambda c : R(m, c)

        self.dimension = m
        self.critical_times = set()
        self.generators = {} # generator_id : time
        self.relations = {} # relation_id : time

        # add generators as zero-simplices
        for key, time in generators.items():
        # for i in range(len(generators)):
            # print(key)
            # print(time)
            # print(Rm(time))
            # print(generators[i])
            self.critical_times.add(Rm(time))
            # self.generators.append(str(i))
            self.generators[str(key)] = Rm(time)

        # add relations as one-simplices
        for key, time in relations.items():
            u = str(key[0])
            v = str(key[1])

            # check valid filtration
            assert(self.generators[u] <= Rm(time) and self.generators[v] <= Rm(time))

            self.relations[(u,v)] = Rm(time)
            id_str = "{}:{}".format(u, v)

            self.critical_times.add(Rm(time))


    def __str__(self):
        pres_str = "M : {};\n".format(self.get_simplicial_complex().simplices())
        critical_times = []
        for time in self.critical_times:
            critical_times.append(str(time))
        critical_str = "C_M := {};\n".format(", ".join(critical_times));
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

    def get_filtration(self, output=False):
        g = self.get_graph();

        for i in self.critical_times:
            print(i)
            vfilt=lambda v : g.vp.time[v] <= i
            efilt=lambda e: g.ep.time[e] <= i
            u = GraphView(g, vfilt=vfilt, efilt=efilt)
            comp, hist = label_components(u)
            print(comp.a)

            # print(comp)

            # if output:
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

m = 1
# Rm = lambda c : R(m, c)
# print(Rm([0]))
# exit()
generators = {
    0 : [0],
    1 : [1],
    2 : [2]
}

# generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
relations = {
    (0, 1) : [4],
    (1, 2) : [5]
}

# m = 2
# generators = strands
# relations = merge_points
# generators = {
#     0 : [0, 0],
#     1 : [1, 1],
#     2 : [2, 2]
# }
#
# # generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
# relations = {
#     (0, 1) : [4, 4],
#     (1, 2) : [5, 5]
# }

presentation = MergePresentation(m, generators, relations)
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
