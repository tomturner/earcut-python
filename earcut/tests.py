import unittest

from earcut.earcut import EarCut


class TestEarCut(unittest.TestCase):

    def test_indices_2d(self):
        e = EarCut()

        indices = e.earcut([10, 0, 0, 50, 60, 60, 70, 10])
        expect_result = [1, 0, 3, 3, 2, 1]

        self.assertEqual(indices, expect_result)

    def test_indices_3d(self):
        e = EarCut()
        indices = e.earcut([10, 0, 0, 0, 50, 0, 60, 60, 0, 70, 10, 0], None, 3)

        expect_result = [1, 0, 3, 3, 2, 1]

        self.assertEqual(indices, expect_result)

    def test_empty(self):
        e = EarCut()
        indices = e.earcut([], [])

        expect_result = []

        self.assertEqual(indices, expect_result)
