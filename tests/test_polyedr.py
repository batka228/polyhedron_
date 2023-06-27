from pytest import approx
from common.r3 import R3
from shadow.polyedr import Segment, Edge, Facet, Polyedr
from math import pi, cos
from functools import reduce
from operator import add
from common.tk_drawer import TkDrawer

class PolyedrEdge:
    
    # При изменении коэф. гомотетии в n раз, 
    # сумма площадей проекций полиэдра изменяется в n^2
    def test_1(self):
        A1 = Polyedr(f"data/box.geom").A()
        A2 = Polyedr(f"data/box1.geom").A()
        assert A1 == approx(A2 / 4)

    # Углы Эйлера влияют на результат
    def test_2(self):
        A1 = Polyedr(f"data/box.geom").A()
        A2 = Polyedr(f"data/box2.geom").A()
        assert A1 != A2

    def test_3(self):
        A = Polyedr(f"data/box.geom").A()
        assert A == approx(75903.68414342678)

    def test_4(self):
        A = Polyedr(f"ccc/.geom").A()
        assert A == approx(0)

    def test_5(self):
        A = Polyedr(f"cube/.geom").A()
        assert A == approx(0)