"""Periodic square-lattice clusters (possibly tilted) for Hubbard ED.

A cluster is defined by two integer superlattice vectors a1, a2. Sites are the
integer points of one fundamental cell. Neighbours and distances use periodic
(minimum-image) boundaries. Conventions follow method/ref/multisize_ED_table_kuku.md:
  2x4   (2,0),(0,4)     sqrt8 (2,2),(-2,2)    sqrt10 (3,1),(-1,3)
  3x4   (3,0),(0,4)     4x4   (4,0),(0,4)
"""
import math
import numpy as np

CLUSTERS = {
    "2x2": ((2, 0), (0, 2)),
    "2x4": ((2, 0), (0, 4)),
    "sqrt8": ((2, 2), (-2, 2)),
    "sqrt10": ((3, 1), (-1, 3)),
    "3x4": ((3, 0), (0, 4)),
    "4x4": ((4, 0), (0, 4)),
}

DIRS = {"+x": (1, 0), "-x": (-1, 0), "+y": (0, 1), "-y": (0, -1)}
# d-wave form factor g(delta): +1 along x, -1 along y (Ku-ku convention)
GFORM = {"+x": 1.0, "-x": 1.0, "+y": -1.0, "-y": -1.0}


class Cluster:
    def __init__(self, name):
        self.name = name
        a1, a2 = CLUSTERS[name]
        self.a1 = np.array(a1, dtype=int)
        self.a2 = np.array(a2, dtype=int)
        self.A = np.array([a1, a2], dtype=float).T  # columns a1, a2
        self.Ainv = np.linalg.inv(self.A)
        self.N = int(round(abs(np.linalg.det(self.A))))
        span = int(np.abs(np.array([a1, a2])).sum()) + 1
        found = {}
        for x in range(-span, span + 1):
            for y in range(-span, span + 1):
                r = self.reduce((x, y))
                found[r] = True
        self.sites = sorted(found)
        assert len(self.sites) == self.N, (name, len(self.sites), self.N)
        self.index = {r: i for i, r in enumerate(self.sites)}
        # neighbours
        self.nbr = {d: [self.index[self.reduce((x + v[0], y + v[1]))] for (x, y) in self.sites]
                    for d, v in DIRS.items()}
        # translations: translation g shifts every site by sites[g]
        self.trans = np.array([[self.index[self.reduce((x + gx, y + gy))] for (x, y) in self.sites]
                               for (gx, gy) in self.sites], dtype=np.int64)

    def reduce(self, r):
        v = np.array(r, dtype=float)
        st = self.Ainv @ v
        st = st - np.floor(st + 1e-9)
        w = self.A @ st
        return (int(round(w[0])), int(round(w[1])))

    def min_image(self, d):
        best = None
        for m in range(-3, 4):
            for n in range(-3, 4):
                v = np.array(d) + m * self.a1 + n * self.a2
                dist = math.hypot(v[0], v[1])
                if best is None or dist < best - 1e-12:
                    best = dist
        return best

    def distance(self, i, j):
        (xi, yi), (xj, yj) = self.sites[i], self.sites[j]
        return self.min_image((xj - xi, yj - yi))

    def r_max(self):
        return max(self.distance(0, j) for j in range(self.N))

    def shell(self, i, R, tol=1e-9):
        return [j for j in range(self.N) if abs(self.distance(i, j) - R) < tol]

    def hop_bonds(self):
        """Directed hops (i -> j) for H = -t sum_i sum_{d in +x,+y} (c_j^dag c_i + h.c.).

        On a length-2 periodic direction, +x and -x lead to the same site, so that
        bond is counted twice (matrix element -2t). This is the plain sum over the
        four neighbour vectors and is checked against Ku-ku's published E0.
        """
        out = []
        for i in range(self.N):
            for d in ("+x", "+y"):
                j = self.nbr[d][i]
                if j == i:
                    continue
                out.append((i, j))
                out.append((j, i))
        return out

    def single_particle_h(self, t=1.0):
        h = np.zeros((self.N, self.N))
        for i, j in self.hop_bonds():
            h[j, i] += -t
        return h

    def k_points(self):
        """Allowed crystal momenta K (2-vectors) with K.a1, K.a2 in 2*pi*Z, reduced to [0, 2pi)."""
        B = 2 * math.pi * np.linalg.inv(self.A).T  # columns b1, b2 with a_i . b_j = 2pi delta
        ks = {}
        for m in range(-self.N, self.N + 1):
            for n in range(-self.N, self.N + 1):
                k = m * B[:, 0] + n * B[:, 1]
                k = np.mod(k, 2 * math.pi)
                key = (round(k[0] / (2 * math.pi) * 1e6) % 1000000, round(k[1] / (2 * math.pi) * 1e6) % 1000000)
                if key not in ks:
                    ks[key] = k
        out = [ks[k] for k in sorted(ks)]
        assert len(out) == self.N, (self.name, len(out))
        return out
