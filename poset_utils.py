
class R():

    def __init__(self, m, coordinate):
        assert(len(coordinate) == m)

        self.dimension = m
        self.coordinate = coordinate

    def copy(self):
        return R(self.dimension, self.coordinate.copy())

    def __hash__(self):
        # return hash(tuple(self.coordinate))
        return hash(tuple(self.coordinate))

    def get_component(self, i):
        assert(1 <= i and i <= self.dimension)
        return self.coordinate[i - 1]

    def to_list(self):
        return self.coordinate

    def shift(self, epsilon):
        m = self.dimension
        e = R(m,  m * [epsilon])
        # print(e)
        return self + e

    def __str__(self):
        components = []
        for i in range(0, self.dimension):
            components.append(str(self.coordinate[i]))
            # text += str(self.coordinate[i])
            # if i < self.dimension - 1:
            #     text += ","
        return "R({})".format(", ".join(components))

    def __add__(self, other):
        assert(self.dimension == other.dimension)
        m = self.dimension
        coordinate = m * [0]
        for i in range(m):
            coordinate[i] = self.coordinate[i] + other.coordinate[i]
        return R(m, coordinate)

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
        for i in range(0, self.dimension):
            if self.coordinate[i] != other.coordinate[i]:
                return False
        return True

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
# print(z >= y)

assert(x.shift(1) == z)
assert(x + x == z)
# print(x)
# print(x.shift(1))
# print(x.shift(1) == z)
# print(x)
# print(y)
# print(x + y)

# assert(x + z == R(2, [3, 3]))

# assert(x.shift(1) == z)

# def downset(collection, basetime):
#     down = []
#     for element in collection:
#         if element <= basetime:
#             down.append(element)
#     return down

def downset(times, basetime):
    down = []
    for i in range(len(times)):
    # for key, time in d.items():
        # print(time)
        # print(basetime)
        time = times[i]
        if time <= basetime:
            down.append(i)
    return down

def upset_times(times, basetime):
    up = []
    for i in range(len(times)):
        time = times[i]
        # print(time)
        # print(basetime)
        # print("time-{}".format(time))
        # print("basetime-{}".format(basetime))
        if time >= basetime:
            # print("--added-time")
            up.append(time)
    return up

# def upset_dict(d, basetime):
#     up = []
#     for key, time in d.items():
#         if time >= basetime:
#             up.append(time)
#     return up

# m = 2
# Rm = lambda c : R(m, c)
# subset = {Rm([0, 0]), Rm([0, 1]), Rm([1, 0])}
# down = downset(subset, Rm([1, 1]))
#
# for element in down:
#     print(element)
