
from poset_utils import *

from itertools import combinations
from simplicial import *
from graph_tool.all import *

DEBUG = False;

class UnionFind:

    def __init__(self, s):
        self.set = s
        self.parent = {}
        for i in s:
            self.parent[i] = i

    def add(self, k):
        self.set.append(k)
        self.parent[k] = k

    def contains(self, k):
        if k in self.set:
            return True
        return False

    def find(self, k):
        if self.parent[k] == k:
            return k
        return self.find(self.parent[k])

    def unite(self, a, b):
        x = self.find(a)
        y = self.find(b)
        self.parent[x] = y

    def equivalent(self, a, b):
        return self.find(a) == self.find(b)

    def __str__(self):
        r = range(len(self.set))
        s = self.set
        return str([[s[j] for j in r if self.equivalent(s[j], s[i])] for i in r if s[i] == self.parent[s[i]]])
        # return str([self.find(i) for i in self.set])

# s = [0, 2, 3, 4, 9, 10]
# print(s)
# u = UnionFind(s)
# u.unite(0, 9)
# print(u)
# u.unite(9, 10)
# print(u)
# print(u.equivalent(0, 10))
# exit()

class MergePoint():
    def __init__(self, id_x, id_y, rm_time):
        self.id_x = id_x
        self.id_y = id_y
        self.rm_time = rm_time

    def __str__(self):
        return "{},{} merge at {}".format(self.id_x, self.id_y, self.rm_time)

class Merge:

    def __init__(self, m, strands, merge_points):
        # return None
        # if DEBUG : print("Merge : init")
        self.dimension = m
        Rm = lambda c : R(m, c)

        self.m_set = {} # M(t) : t in C_M.

        # generate critical times set
        self.critical_times = set()
        for index, time in strands.items():
            self.critical_times.add(Rm(time))
        for indices, times in merge_points.items():
            for time in times:
                self.critical_times.add(Rm(time))

        # initialize internal strands
        self.strand_times = []
        self.strand_ids = {}
        for key, time in strands.items():
            index = len(self.strand_times)
            self.strand_times.append(Rm(time))
            self.strand_ids[key] = index

        # initialize sets
        for time in self.critical_times:
            self.m_set[time] = UnionFind(downset(self.strand_times, time))

        # initialize internal merge points tracker
        self.merge_points = []
        for key, times in merge_points.items():
            for i in range(0, len(key) - 1):
                id_x = self.strand_ids[key[i]]
                id_y = self.strand_ids[key[i + 1]]
                for time in times:
                    m_point = MergePoint(id_x, id_y, Rm(time))
                    self.merge_points.append(m_point)
            # print("{} :: {}".format(key, times))
            # print(time)
        for m_point in self.merge_points:
            print("--{}".format(m_point))

        # merge sets
        for key, times in merge_points.items():
            for time in times:
                uptimes = upset_times(list(self.critical_times), Rm(time))
                for u_time in uptimes:
                    for i in range(0, len(key) - 1):
                        current_id = self.strand_ids[key[i]]
                        next_id = self.strand_ids[key[i + 1]]

                        self.m_set[u_time].unite(current_id, next_id)

    def add_critical_time(self, time):
        self.critical_times.add(R(self.dimension, time))

    def add_strand_time(self, time):
        self.strand_times.append(R(self.dimension, time))

    def add_strand(self, id, birth_time):
        # check id is unique
        assert(id not in self.strand_ids)
        rm_time = R(self.dimension, birth_time)


        # update internal strands
        strand_index = len(self.strand_times)
        self.strand_ids[id] = strand_index

        self.add_strand_time(birth_time)
        self.add_critical_time(birth_time)

        # update m_sets[t] where t : upset(time)
        if rm_time not in self.m_set:
            self.add_m_set(birth_time)
            # self.m_set[rm_time] = UnionFind(downset(self.strand_times, rm_time))

        for time in self.critical_times:
            if time > rm_time:
                self.m_set[time].add(strand_index)

        return strand_index

    def add_m_set(self, time):
        # update m_set[time] to have all strands alive at `time`
        rm_time = R(self.dimension, time)
        self.m_set[rm_time] = UnionFind(downset(self.strand_times, rm_time))

        # merge points identified before `time`
        for m_point in self.merge_points:
            if m_point.rm_time <= rm_time:
                self.m_set[rm_time].unite(m_point.id_x, m_point.id_y)


    # def merge(self, merge_time):
    #     # Merge everything at merge_time
    #     for merge_point in self.merge_points:
    #         print(merge_point)
    #
    #     return None

    def add_merge_point(self, key_x, key_y, merge_time):
        # m = self.dimension
        rm_time = R(self.dimension, merge_time)


        id_x = self.strand_ids[key_x]
        id_y = self.strand_ids[key_y]
        # check both born <= merge_time
        assert(self.strand_times[id_x] <= rm_time)
        assert(self.strand_times[id_y] <= rm_time)

        # add merge_time as critical time if doesn't exist
        if rm_time not in self.m_set:
            self.add_critical_time(merge_time)
            # self.critical_times.add(R(m, merge_time))
            self.add_m_set(merge_time)
            # self.merge(merge_time)

        # add merge point to internal tracker
        m_point = MergePoint(id_x, id_y, rm_time)
        self.merge_points.append(m_point)


        uptimes = upset_times(list(self.critical_times), rm_time)
        for u_time in uptimes:
            self.m_set[u_time].unite(id_x, id_y)

        return None

    def count_strands(self):
        return len(self.strand_times)

    def get_id_by_index(self, index):
        for key, value in self.strand_ids.items():
            if value == index:
                return key
        return None

    def get_strand_by_index(self, index):
        # key = self.strand_indices[index]
        # birth_time = self.strands[key]
        # return None
        # key = self.strand_ids[index]
        key = self.get_id_by_index(index)
        birth_time = self.strand_times[index]
        return "--{} : strands[{}] = {}".format(index, key, birth_time)

    def get_set(self, r):
        return None

    def are_equal(self, s, t, i, j):
        """
        Returns if M(s <= t)(g_i) == M(s <= t)(g_j)
        """
        assert(s <= t)
        assert(i in self.m_set[s].set)
        assert(j in self.m_set[s].set)
        return self.m_set[t].equivalent(i, j)

    def get_structure_map(self, s, t): # M(s) to M(t) inclusion
        assert(s <= t)
        return None

    # def evaluate_structure_map(self, s, t, i):
    #     assert(s <= t)
    #     assert(i in self.m_set[s].set)
    #     # print(self.m_set[t].equivalent(0, 1))
    #     return None

    def shift(self, epsilon):
        epsilon = Rm(self.dimension * [epsilon])

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
        components = []
        for key, value in self.m_set.items():
            components.append("M({}) = {}".format(key, value))
        m_sets = "{}".format("\n".join(components))
            # print()

        # build critical times string
        components = []
        for time in self.critical_times:
            components.append(str(time))
        critical_times = "C_M = [{}]".format(", ".join(components)) + "\n"

        text = critical_times + m_sets + "\n"
        return text


m = 2
Rm = lambda c : R(m, c)
strands = {
    0 : [0, 0],
    1 : [0, 0],
    22 : [1, 1]
}

# generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
merge_points = {
    (0, 1) : [[0, 1], [1, 0]],
    (0, 1, 22) : [[2, 2]]
}

# merge = Merge(m, strands, merge_points)
# strand_count = merge.count_strands()
# for i in range(strand_count):
#     print(merge.get_strand_by_index(i))
# poset = merge.get_structure_poset()
# for pair in poset:
#     print("{} <= {}".format(pair[0], pair[1]))
    # print(pair[0])
    # print(pair[1])
# ev = merge.evaluate_structure_map(Rm([0, 0]), Rm([0, 1]), 0)
# print(ev)
# print(merge.are_equal(Rm([0, 0]), Rm([0, 1]), 0, 1))
# print(merge)

m = 1
Rm = lambda c : R(m, c)

strands = {
    0 : [0],
    1 : [1],
    2 : [2]
}
merge_points = {
    (0, 1, 2) : [[3]]
}
M = Merge(m, strands, merge_points)
# assert(M.count_strands() == 3)
print(M)
M.add_strand("3", [1.5])
print(M)
strand_count = M.count_strands()
for i in range(strand_count):
    print(M.get_strand_by_index(i))# M.add_merge_point
M.add_merge_point(0, "3", [5])
print(M)

strands = {0 : [2]}
merge_points = {}
N = Merge(m, strands, merge_points)
# assert(N.count_strands() == 1)
print(N)

class MergeMorphism:

    def __init__(self, M, N, phi):
        assert(M.dimension == N.dimension)

        self.source = M
        self.target = N
        self.morphism = phi

    def is_valid(self):
        return False

phi = MergeMorphism(M, N, [0, 0, 0])



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

    def is_minimal(self):

        return False

    def get_critical_times(self):
        return self.critical_times;

    def get_critical_times_i(self, i):
        # assert(1 <= i and i <= self.dimension)
        critical_times_i = set()
        for element in self.critical_times:
            critical_times_i.add(element.get_component(i))
        return critical_times_i

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

def get_possible_interleaving_values(M, N):
    assert(M.dimension == N.dimension)

    m = M.dimension

    cm = []
    cn = []

    for i in range(1, m + 1):
        cm.append(M.get_critical_times_i(i))
        cn.append(N.get_critical_times_i(i))

    possible_distances = set()
    for i in range(0, m):
        # TODO : rewrite this to use `combinations`
        for x in cm[i]:
            for y in cn[i]:
                possible_distances.add(abs(x - y))
        for x, y in combinations(cm[i], 2):
            possible_distances.add(0.5 * abs(x - y))
        for x, y in combinations(cn[i], 2):
            possible_distances.add(0.5 * abs(x - y))

    return possible_distances

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
    (0, 1) : [3],
    (1, 2) : [4]
}

PM = MergePresentation(m, generators, relations)

generators = {0 : [2]}
relations = {}

PN = MergePresentation(m, generators, relations)

d = get_possible_interleaving_values(PM, PN)
# print(d)
#
# exit()

m = 2

# generators = strands
# relations = merge_points
generators = {
    0 : [0, 0],
    1 : [1, 1],
    2 : [2, 2]
}

# generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
relations = {
    (0, 1) : [4, 4],
    (1, 2) : [5, 5]
}
presentationA = MergePresentation(m, generators, relations)


# generators = strands
# relations = merge_points
generators = {
    0 : [0, 1],
    1 : [1, 2],
    2 : [2, 3]
}

# generators = [0, 1, 2] # g_0 born at 0, g_1 born at 0, g_2 born at 0
relations = {
    (0, 1) : [4, 7],
    (1, 2) : [5, 6]
}

presentationB = MergePresentation(m, generators, relations)
# print(presentation)
# print(presentation.count_connected_components())
# print(presentation.is_connected())
# g = presentation.get_graph()
# presentation.get_filtration()
# critical_times_i = presentation.get_critical_times_i(2)
# print(critical_times_i)

d = get_possible_interleaving_values(presentationA, presentationB)
print(d)

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
