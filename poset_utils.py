
class R():

    def __init__(self, m, coordinate):
        assert(len(coordinate) == m)

        self.dimension = m
        self.coordinate = coordinate

    def __hash__(self):
        # return hash(tuple(self.coordinate))
        return hash(tuple(self.coordinate))

    def __str__(self):
        components = []
        for i in range(0, self.dimension):
            components.append(str(self.coordinate[i]))
            # text += str(self.coordinate[i])
            # if i < self.dimension - 1:
            #     text += ","
        return "R({})".format(", ".join(components))

    def __lt__(self, other):
        less_than = True
        exist_strict = False
        for i in range(0, self.dimension):
            # print((self.coordinate[i] <= other.coordinate[i]))
            less_than = less_than and (self.coordinate[i] <= other.coordinate[i])
            if self.coordinate[i] < other.coordinate[i]:
                exist_strict = True

        return less_than and exist_strict

    def __eq__(self, other):
        equal = True
        for i in range(0, self.dimension):
            if self.coordinate[i] != other.coordinate[i]:
                equal = False
        return equal

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return self > other or self == other

# Unit Tests

x = R(2, [1, 1])
y = R(2, [1, 2])
z = R(2, [2, 2])

assert(x < z)
assert(x <= z)
assert(z > x)
assert(z >= x)
assert(x < y)
assert(x <= y)
assert(z > y)
assert(z >= y)
assert(x <= x)
assert(not x < x)
assert(not x > x)

def downset(collection, basetime):
    down = set()
    for element in collection:
        if element <= basetime:
            down.add(element)
    return down

def downset_dict(d, basetime):
    down = set()
    for key, time in d.items():
        # print(time)
        # print(basetime)
        if time <= basetime:
            down.add(key)
    return down

# m = 2
# Rm = lambda c : R(m, c)
# subset = {Rm([0, 0]), Rm([0, 1]), Rm([1, 0])}
# down = downset(subset, Rm([1, 1]))
#
# for element in down:
#     print(element)
