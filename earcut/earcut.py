import math

# look at https://github.com/tomturner/earcut-python
# mirrors JS version https://github.com/mapbox/earcut
# version 2.1.5


class Node(object):

    __slots__ = ('i',
                 'x',
                 'y',
                 'prev',
                 'next',
                 'z',
                 'prev_z',
                 'next_z',
                 'steiner')

    def __init__(self, i, x, y):
        # vertex index in coordinates array
        self.i = i

        # vertex coordinates

        self.x = x
        self.y = y

        # previous and next vertex nodes in a polygon ring
        self.prev = None
        self.next = None

        # z-order curve value
        self.z = None

        # previous and next nodes in z-order
        self.prev_z = None
        self.next_z = None

        # indicates whether this is a steiner point
        self.steiner = False


class EarCut(object):
    def earcut(self, data, hole_indices=None, dim=None):
        dim = dim or 2

        has_holes = hole_indices and len(hole_indices)
        outer_len = hole_indices[0] * dim if has_holes else len(data)
        outer_node = self.linked_list(data, 0, outer_len, dim, True)
        triangles = []

        if not outer_node or outer_node.next == outer_node.prev:
            return triangles

        min_x = None
        min_y = None
        inv_size = None

        if has_holes:
            outer_node = self.eliminate_holes(data, hole_indices, outer_node, dim)

        # if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
        if len(data) > 80 * dim:
            min_x = max_x = data[0]
            min_y = max_y = data[1]

            for i in range(dim, outer_len, dim):
                x = data[i]
                y = data[i + 1]
                if x < min_x:
                    min_x = x
                if y < min_y:
                    min_y = y
                if x > max_x:
                    max_x = x
                if y > max_y:
                    max_y = y

            # min_x, min_y and size are later used to transform coords into integers for z-order calculation
            inv_size = max(max_x - min_x, max_y - min_y)

            inv_size = 1 / inv_size if inv_size != 0 else 0

        self.earcut_linked(outer_node, triangles, dim, min_x, min_y, inv_size)

        return triangles

    def linked_list(self, data, start, end, dim, clockwise):
        """
        create a circular doubly linked _list from polygon points in the specified winding order
        """
        last = None

        if clockwise == (self.signed_area(data, start, end, dim) > 0):
            for i in range(start, end, dim):
                last = self.insert_node(i, data[i], data[i + 1], last)

        else:
            for i in reversed(range(start, end, dim)):
                last = self.insert_node(i, data[i], data[i + 1], last)

        if last and self.equals(last, last.next):
            self.remove_node(last)
            last = last.next

        return last

    def filter_points(self, start, end=None):
        """
        eliminate colinear or duplicate points
        """
        if not start:
            return start
        if not end:
            end = start

        p = start
        again = True

        while again or p != end:
            again = False

            if not p.steiner and (self.equals(p, p.next) or self.area(p.prev, p, p.next) == 0):
                self.remove_node(p)
                p = end = p.prev
                if p == p.next:
                    break

                again = True

            else:
                p = p.next

        return end

    def earcut_linked(self, ear, triangles, dim, min_x, min_y, inv_size, _pass=None):
        """
        main ear slicing loop which triangulates a polygon (given as a linked _list)
        """
        if not ear:
            return

        # interlink polygon nodes in z-order
        if not _pass and inv_size:
            self.index_curve(ear, min_x, min_y, inv_size)

        stop = ear

        # iterate through ears, slicing them one by one
        while ear.prev != ear.next:
            _prev = ear.prev
            _next = ear.next

            if self.is_ear_hashed(ear, min_x, min_y, inv_size) if inv_size else self.is_ear(ear):
                # cut off the triangle
                triangles.append(_prev.i // dim)
                triangles.append(ear.i // dim)
                triangles.append(_next.i // dim)

                self.remove_node(ear)

                # skipping the next vertex leads to less sliver triangles
                ear = _next.next
                stop = _next.next

                continue

            ear = _next

            # if we looped through the whole remaining polygon and can't find any more ears
            if ear == stop:
                # try filtering points and slicing again
                if not _pass:
                    self.earcut_linked(self.filter_points(ear), triangles, dim, min_x, min_y, inv_size, 1)

                    # if this didn't work, try curing all small self-intersections locally
                elif _pass == 1:
                    ear = self.cure_local_intersections(ear, triangles, dim)
                    self.earcut_linked(ear, triangles, dim, min_x, min_y, inv_size, 2)

                    # as a last resort, try splitting the remaining polygon into two
                elif _pass == 2:
                    self.split_earcut(ear, triangles, dim, min_x, min_y, inv_size)

                break

    def is_ear(self, ear):
        """
        check whether a polygon node forms a valid ear with adjacent nodes
        """
        a = ear.prev
        b = ear
        c = ear.next

        if self.area(a, b, c) >= 0:
            return False  # reflex, can't be an ear

        # now make sure we don't have other points inside the potential ear
        p = ear.next.next

        while p != ear.prev:
            if self.point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) and self.area(p.prev, p, p.next) >= 0:
                return False
            p = p.next

        return True

    def is_ear_hashed(self, ear, min_x, min_y, inv_size):
        a = ear.prev
        b = ear
        c = ear.next

        if self.area(a, b, c) >= 0:
            return False  # reflex, can't be an ear

        # triangle bbox; min & max are calculated like this for speed
        min_t_x = (a.x if a.x < c.x else c.x) if a.x < b.x else (b.x if b.x < c.x else c.x)
        min_t_y = (a.y if a.y < c.y else c.y) if a.y < b.y else (b.y if b.y < c.y else c.y)
        max_t_x = (a.x if a.x > c.x else c.x) if a.x > b.x else (b.x if b.x > c.x else c.x)
        max_t_y = (a.y if a.y > c.y else c.y) if a.y > b.y else (b.y if b.y > c.y else c.y)

        # z-order range for the current triangle bbox;
        min_z = self.z_order(min_t_x, min_t_y, min_x, min_y, inv_size)
        max_z = self.z_order(max_t_x, max_t_y, min_x, min_y, inv_size)

        p = ear.prev_z
        n = ear.next_z

        # look for points inside the triangle in both directions
        while p and p.z >= min_z and n and n.z <= max_z:
            if p != ear.prev and p != ear.next and self.point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x,
                                                                          p.y) and self.area(p.prev, p, p.next) >= 0:
                return False
            p = p.next_z
            if (n != ear.prev and n != ear.next and self.point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) and
                    self.area(n.prev, n, n.next) >= 0):
                return False
            n = n.nextZ

        # look for remaining points in decreasing z-order
        while p and p.z >= min_z:
            if p != ear.prev and p != ear.next and self.point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x,
                                                                          p.y) and self.area(p.prev, p, p.next) >= 0:
                return False
            p = p.prev_z

        # look for remaining points in increasing z-order
        while n and n.z <= max_z:
            if (n != ear.prev and n != ear.next and
                self.point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) and
                    self.area(n.prev, n, n.next) >= 0):
                return False
            n = n.next_z

        return True

    def cure_local_intersections(self, start, triangles, dim):
        """
        go through all polygon nodes and cure small local self-intersections
        """
        do = True
        p = start

        while do or p != start:
            do = False

            a = p.prev
            b = p.next.next

            if (not self.equals(a, b) and
                    self.intersects(a, p, p.next, b) and
                    self.locally_inside(a, b) and
                    self.locally_inside(b, a)):
                triangles.append(a.i // dim)
                triangles.append(p.i // dim)
                triangles.append(b.i // dim)

                # remove two nodes involved
                self.remove_node(p)
                self.remove_node(p.next)

                p = start = b

            p = p.next

        return p

    def split_earcut(self, start, triangles, dim, min_x, min_y, inv_size):
        """
        try splitting polygon into two and triangulate them independently
        look for a valid diagonal that divides the polygon into two
        """
        do = True
        a = start

        while do or a != start:
            do = False
            b = a.next.next

            while b != a.prev:
                if a.i != b.i and self.is_valid_diagonal(a, b):
                    # split the polygon in two by the diagonal
                    c = self.split_polygon(a, b)

                    # filter colinear points around the cuts
                    a = self.filter_points(a, a.next)
                    c = self.filter_points(c, c.next)

                    # run earcut on each half
                    self.earcut_linked(a, triangles, dim, min_x, min_y, inv_size)
                    self.earcut_linked(c, triangles, dim, min_x, min_y, inv_size)
                    return

                b = b.next

            a = a.next

    def eliminate_holes(self, data, hole_indices, outer_node, dim):
        """
        link every hole into the outer loop, producing a single-ring polygon without holes
        """
        queue = []
        _len = len(hole_indices)

        for i in range(len(hole_indices)):
            start = hole_indices[i] * dim
            end = hole_indices[i + 1] * dim if i < _len - 1 else len(data)
            _list = self.linked_list(data, start, end, dim, False)

            if _list == _list.next:
                _list.steiner = True

            queue.append(self.get_leftmost(_list))

        queue = sorted(queue, key=lambda _i: _i.x)

        # process holes from left to right
        for i in range(len(queue)):
            self.eliminate_hole(queue[i], outer_node)
            outer_node = self.filter_points(outer_node, outer_node.next)

        return outer_node

    @staticmethod
    def compare_x(a, b):
        return a.x - b.x

    def eliminate_hole(self, hole, outer_node):
        """
        find a bridge between vertices that connects hole with an outer ring and and link it
        """
        outer_node = self.find_hole_bridge(hole, outer_node)
        if outer_node:
            b = self.split_polygon(outer_node, hole)
            self.filter_points(b, b.next)

    def find_hole_bridge(self, hole, outer_node):
        """§
        David Eberly's algorithm for finding a bridge between hole and outer polygon
        """
        do = True
        p = outer_node
        hx = hole.x
        hy = hole.y
        qx = -math.inf
        m = None

        # find a segment intersected by a ray from the hole's leftmost point to the left;
        # segment's endpoint with lesser x will be potential connection point
        while do or p != outer_node:
            do = False
            if hy <= p.y and hy >= p.next.y != p.y:
                x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y)

                if hx >= x > qx:
                    qx = x

                    if x == hx:
                        if hy == p.y:
                            return p
                        if hy == p.next.y:
                            return p.next

                    m = p if p.x < p.next.x else p.next

            p = p.next

        if not m:
            return None

        if hx == qx:
            return m.prev  # hole touches outer segment; pick lower endpoint

        # look for points inside the triangle of hole point, segment intersection and endpoint;
        # if there are no points found, we have a valid connection;
        # otherwise choose the point of the minimum angle with the ray as connection point

        stop = m
        mx = m.x
        my = m.y
        tan_min = math.inf

        p = m.next

        while p != stop:
            hx_or_qx = hx if hy < my else qx
            qx_or_hx = qx if hy < my else hx

            if hx >= p.x >= mx and hx != p.x and self.point_in_triangle(hx_or_qx, hy, mx, my, qx_or_hx, hy, p.x, p.y):

                tan = abs(hy - p.y) / (hx - p.x)  # tangential

                if (tan < tan_min or (tan == tan_min and p.x > m.x)) and self.locally_inside(p, hole):
                    m = p
                    tan_min = tan

            p = p.next

        return m

    def index_curve(self, start, min_x, min_y, inv_size):
        """
        interlink polygon nodes in z-order
        """
        do = True
        p = start

        while do or p != start:
            do = False

            if p.z is None:
                p.z = self.z_order(p.x, p.y, min_x, min_y, inv_size)

            p.prev_z = p.prev
            p.prev_z = p.next
            p = p.next

        p.prev_z.next_z = None
        p.prev_z = None

        self.sort_linked(p)

    @staticmethod
    def sort_linked(_list):
        """
        Simon Tatham's linked _list merge sort algorithm
        http:#www.chiark.greenend.org.uk/~sgtatham/algorithms/_listsort.html
        """
        do = True
        num_merges = None
        in_size = 1

        while do or num_merges > 1:
            do = False
            p = _list
            _list = None
            tail = None
            num_merges = 0

            while p:
                num_merges += 1
                q = p
                p_size = 0
                for i in range(in_size):
                    p_size += 1
                    q = q.next_z
                    if not q:
                        break

                q_size = in_size

                while p_size > 0 or (q_size > 0 and q):

                    if p_size != 0 and (q_size == 0 or not q or p.z <= q.z):
                        e = p
                        p = p.next_z
                        p_size -= 1

                    else:
                        e = q
                        q = q.next_z
                        q_size -= 1

                    if tail:
                        tail.next_z = e

                    else:
                        _list = e

                    e.prev_z = tail
                    tail = e

                p = q

            tail.next_z = None
            in_size *= 2

        return _list

    @staticmethod
    def z_order(x, y, min_x, min_y, inv_size):
        """
        z-order of a point given coords and inverse of the longer side of data bbox
        coords are transformed into non-negative 15-bit integer range
        """
        #
        x = int(32767 * (x - min_x) * inv_size)
        y = int(32767 * (y - min_y) * inv_size)

        x = (x | (x << 8)) & 0x00FF00FF
        x = (x | (x << 4)) & 0x0F0F0F0F
        x = (x | (x << 2)) & 0x33333333
        x = (x | (x << 1)) & 0x55555555

        y = (y | (y << 8)) & 0x00FF00FF
        y = (y | (y << 4)) & 0x0F0F0F0F
        y = (y | (y << 2)) & 0x33333333
        y = (y | (y << 1)) & 0x55555555

        return x | (y << 1)

    @staticmethod
    def get_leftmost(start):
        """
        find the leftmost node of a polygon ring
        """
        do = True
        p = start
        leftmost = start

        while do or p != start:
            do = False
            if p.x < leftmost.x or (p.x == leftmost.x and p.y < leftmost.y):
                leftmost = p
            p = p.next

        return leftmost

    @staticmethod
    def point_in_triangle(ax, ay, bx, by, cx, cy, px, py):
        """
        check if a point lies within a convex triangle
        """
        return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 and \
               (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 and \
               (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0

    def is_valid_diagonal(self, a, b):
        """
        check if a diagonal between two polygon nodes is valid (lies in polygon interior)
        """
        return (a.next.i != b.i and a.prev.i != b.i and not self.intersects_polygon(a, b) and
                self.locally_inside(a, b) and self.locally_inside(b, a) and self.middle_inside(a, b))

    @staticmethod
    def area(p, q, r):
        """
        signed area of a triangle
        """
        return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

    @staticmethod
    def equals(p1, p2):
        """
        check if two points are equal
        """
        return p1.x == p2.x and p1.y == p2.y

    def intersects(self, p1, q1, p2, q2):
        """
        check if two segments intersect
        """
        if (self.equals(p1, q1) and self.equals(p2, q2)) or (self.equals(p1, q2) and self.equals(p2, q1)):
            return True

        return (self.area(p1, q1, p2) > 0 != self.area(p1, q1, q2) > 0 and
                self.area(p2, q2, p1) > 0 != self.area(p2, q2, q1) > 0)

    def intersects_polygon(self, a, b):
        """
        check if a polygon diagonal intersects any polygon segments
        """
        do = True
        p = a

        while do or p != a:
            do = False
            if (p.i != a.i and p.next.i != a.i and p.i != b.i and
                    p.next.i != b.i and self.intersects(p, p.next, a, b)):
                return True

            p = p.next

        return False

    def locally_inside(self, a, b):
        """
        check if a polygon diagonal is locally inside the polygon
        """
        if self.area(a.prev, a, a.next) < 0:
            return self.area(a, b, a.next) >= 0 and self.area(a, a.prev, b) >= 0
        else:
            return self.area(a, b, a.prev) < 0 or self.area(a, a.next, b) < 0

    @staticmethod
    def middle_inside(a, b):
        """
        check if the middle point of a polygon diagonal is inside the polygon
        """
        do = True
        p = a
        inside = False
        px = (a.x + b.x) / 2
        py = (a.y + b.y) / 2

        while do or p != a:
            do = False
            if (((p.y > py) != (p.next.y > py)) and p.next.y != p.y and
                    (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x)):
                inside = not inside

            p = p.next

        return inside

    @staticmethod
    def split_polygon(a, b):
        """
        link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
        if one belongs to the outer ring and another to a hole, it merges it into a single ring
        """
        a2 = Node(a.i, a.x, a.y)
        b2 = Node(b.i, b.x, b.y)
        an = a.next
        bp = b.prev

        a.next = b
        b.prev = a

        a2.next = an
        an.prev = a2

        b2.next = a2
        a2.prev = b2

        bp.next = b2
        b2.prev = bp

        return b2

    @staticmethod
    def insert_node(i, x, y, last):
        """
        create a node and optionally link it with previous one (in a circular doubly linked _list)
        """
        p = Node(i, x, y)

        if not last:
            p.prev = p
            p.next = p

        else:
            p.next = last.next
            p.prev = last
            last.next.prev = p
            last.next = p

        return p

    @staticmethod
    def remove_node(p):
        p.next.prev = p.prev
        p.prev.next = p.next

        if p.prev_z:
            p.prev_z.next_z = p.next_z

        if p.next_z:
            p.next_z.prev_z = p.prev_z

    def deviation(self, data, hole_indices, dim, triangles):
        """
        return a percentage difference between the polygon area and its triangulation area;
        used to verify correctness of triangulation
        """
        _len = len(hole_indices)
        has_holes = hole_indices and len(hole_indices)
        outer_len = hole_indices[0] * dim if has_holes else len(data)

        polygon_area = abs(self.signed_area(data, 0, outer_len, dim))

        if has_holes:
            for i in range(_len):
                start = hole_indices[i] * dim
                end = hole_indices[i + 1] * dim if i < _len - 1 else len(data)
                polygon_area -= abs(self.signed_area(data, start, end, dim))

        triangles_area = 0

        for i in range(0, len(triangles), 3):
            a = triangles[i] * dim
            b = triangles[i + 1] * dim
            c = triangles[i + 2] * dim
            triangles_area += abs(
                (data[a] - data[c]) * (data[b + 1] - data[a + 1]) -
                (data[a] - data[b]) * (data[c + 1] - data[a + 1]))

        if polygon_area == 0 and triangles_area == 0:
            return 0

        return abs((triangles_area - polygon_area) / polygon_area)

    @staticmethod
    def signed_area(data, start, end, dim):
        _sum = 0
        j = end - dim

        for i in range(start, end, dim):
            _sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1])
            j = i

        return _sum

    # turn a polygon in a multi-dimensional array form (e.g. as in GeoJSON) into a form Earcut accepts
    @staticmethod
    def flatten(data):
        dim = len(data[0][0])
        result = {
            'vertices': [],
            'holes': [],
            'dimensions': dim
        }
        hole_index = 0

        for i in range(len(data)):
            for j in range(len(data[i])):
                for d in range(dim):
                    result['vertices'].append(data[i][j][d])

            if i > 0:
                hole_index += len(data[i - 1])
                result['holes'].append(hole_index)

        return result

    @staticmethod
    def unflatten(data):
        result = []

        for i in range(0, len(data), 3):
            result.append(tuple(data[i:i + 3]))

        return result
