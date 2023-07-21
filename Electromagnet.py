import numpy as np
import scipy.constants
import sympy
from sympy import assuming, Q, N
from sympy.plotting.plot import *
from sympy import sin, cos, exp, pi, symbols
from sympy.vector import CoordSys3D, ParametricRegion, ImplicitRegion, vector_integrate, divergence, curl
from sympy.abc import r, theta, phi, s, t, u

from sympy.utilities.autowrap import ufuncify

sympy.init_printing(use_latex=True)

x = sympy.Symbol('x', real=True, imaginary=False)
y = sympy.Symbol('y', real=True, imaginary=False)
z = sympy.Symbol('z', real=True, imaginary=False)
C = CoordSys3D('C')

class System:
    def __init__(self):
        self.entities = []

    def permittivity_at_point(self, x,y,z):
        return 8.85e-12
        #return scipy.constants.epsilon_0
        #return 1

    def permeability_at_point(self, x, y, z):
        return scipy.constants.mu_0

    def add_entity(self, entity):
        assert isinstance(entity, Entity)
        self.entities.append(entity)

    def calculate_total_field_at_point(self, x, y, z):
        E_total = np.zeros_like([x, y, z], dtype=float)
        for ent in self.entities:
            E = np.nan_to_num(ent.field_at_point(x, y, z))
            E_total += E
        return E_total

    def calculate_total_mag_field_at_point(self, x, y, z):
        B_total = np.zeros_like([x,y,z], dtype=float)
        for ent in self.entities:
            B = np.nan_to_num(ent.magnetic_field_at_point(x, y, z))
            B_total += B
        return B_total

    def plot_entities(self, mpl_ax):
        for ent in self.entities:
            ent.plot_surface(mpl_ax)

    def plot_electric_field(self, mpl_ax):
        #vfunc = np.vectorize(self.calculate_total_field_at_point)
        #return vfunc(X, Y, Z)
        #X,Y,Z = np.arange(10), np.arange(10), np.arange(10)
        X, Y, Z = np.meshgrid(np.arange(-1, 1, 0.4),
                      np.arange(-1, 1, 0.4),
                      np.arange(-1, 1, 0.4))
        E_field = self.calculate_total_field_at_point(X,Y,Z)
        mpl_ax.quiver(X,Y,Z, *E_field, normalize=True, length=0.2)

    def plot_magnetic_field(self, mpl_ax):
        X, Y, Z = np.meshgrid(np.arange(-1, 1, 0.4),
                      np.arange(-1, 1, 0.4),
                      np.arange(-1, 1, 0.4))
        B_field = self.calculate_total_mag_field_at_point(X, Y, Z)
        mpl_ax.quiver(X,Y,Z, *B_field, normalize=True, length=0.2, color='g')

class Entity:
    def __init__(self, system, region):
        self.region = region
        self.system = system

    def field_at_point(self, i,j,k):
        assert False, "This class sucks"
        P = i*C.i + j*C.j + k*C.k
        X = x*C.i + y*C.j + z*C.k
        print(sympy.sqrt(dx**2).is_positive)
        norm = sympy.Abs((P-X).magnitude())
        print(norm)
        R = (P-X)/((P-X).magnitude()**3)
        #R /= 4*sympy.pi*self.system.permittivity_at_point(x, y, z)
        #print(P)
        #print(R)
        #print(self.region)
        return vector_integrate(R, self.region)

class Point(Entity):
    def __init__(self, system, pos, charge):
        Entity.__init__(self, system, None)
        self.pos = pos
        self.charge = charge

    def field_at_point(self, px, py, pz):
        #if self.pos == (px, py, pz): return np.array([0,0,0])
        i,j,k = sympy.symbols("i j k")
        P = i*C.i + j*C.j + k*C.k
        X = self.pos[0]*C.i + self.pos[1]*C.j + self.pos[2]*C.k
        R = self.charge * (P-X)/((P-X).magnitude()**3)
        R /= (4*sympy.pi*self.system.permittivity_at_point(i,j,k))

        E = ufuncify([i,j,k], list(R.to_matrix(C)))
        return np.array(E(px, py, pz))
        #return np.array([N(R.coeff(co), 3) for co in [C.i, C.j, C.k]])

    def plot_surface(self, ax):
        ax.plot(*self.pos, marker="o")

    def force_experienced(self):
        E = self.system.calculate_total_field_at_point(*self.pos)
        return np.multiply(E, self.charge)

class Line(Entity):
    # fx, fy, fz are all functions of t only
    def __init__(self, system, fx, fy, fz, charge_density=1, current=0, t_limits=(0,1)):
        self.system = system
        self.X = 1*fx*C.i + 1*fy*C.j + 1*fz*C.k
        self.charge_density=1
        self.current = current
        self.t_limits = t_limits

    # Electric field at px, py, pz, takes np array
    def field_at_point(self, px, py, pz):
        i,j,k = symbols("i j k")
        P = 1*i*C.i + 1*j*C.j + 1*k*C.k
        R = self.charge_density * (P-self.X)/((P-self.X).magnitude()**3)
        R /= (4*sympy.pi * self.system.permittivity_at_point(0,0,0))

        E = []
        for coeff in list(R.to_matrix(C)):
            dE = ufuncify([t, i,j,k], coeff)
            min_t, max_t = self.t_limits
            E_integral = lambda a,b,c: scipy.integrate.quad(dE, min_t,max_t, args=(a,b,c,), epsrel=0.0001)[0]
            E.append(np.array(np.vectorize(E_integral)(px, py, pz)))
        return E

    # Magnetic field at px, py, pz, also takes np arrays
    def magnetic_field_at_point(self, px, py, pz):
        i,j,k = symbols("i j k")
        P = 1*i*C.i + 1*j*C.j + 1*k*C.k
        # biot-savart law
        R = self.current * self.system.permeability_at_point(px, py, pz)/(4*pi)
        R *= self.X.cross(P)/((P-self.X).magnitude()**2)

        B = []
        for coeff in list(R.to_matrix(C)):
            dB = ufuncify([t, i,j,k], coeff)
            min_t, max_t = self.t_limits
            B_integral = lambda a,b,c: scipy.integrate.quad(dB, min_t,max_t, args=(a,b,c,), epsrel=0.0001)[0]
            B.append(np.array(np.vectorize(B_integral)(px, py, pz)))
        return B

    # Pass a matplotlib axis to plot self
    def plot_surface(self, ax):
        t = sympy.symbols('t')
        line = sympy.lambdify([t], self.X.to_matrix(C))
        min_t, max_t = self.t_limits
        min_t = np.nan_to_num(float(min_t), posinf=10e5, neginf=-10e5)
        max_t = np.nan_to_num(float(max_t), posinf=10e5, neginf=-10e5)

        X = []
        Y = []
        Z = []
        for dt in np.linspace(min_t, max_t, 20):
            o = line(dt)
            X.append(o[0])
            Y.append(o[1])
            Z.append(o[2])
        ax.plot(X,Y,Z, color='orange')

# Parametric coordiantes each are functions of surfaces of form (t, u)
# charge_density must also be an sympy expression of (t, u)
class Surface(Entity):
    def __init__(self, system, fx, fy, fz, charge_density=1):
        self.system = System
        self.X = fx*C.i + fy*C.j + fz*C.k
        #print(type(self.X))
        self.charge_density=charge_density

    def field_at_point(self, px, py, pz):
        i,j,k = sympy.symbols("i j k", real=True)

        P = i*C.i + j*C.j + k*C.k # point subbed in
        R = (P-self.X)/(P-self.X).magnitude()**3 # directional com
        R *= self.charge_density / (4*sympy.pi*system.permittivity_at_point(0,0,0))
        E = []
        for coeff in list(R.to_matrix(C)):
            dE = ufuncify([u,t,i,j,k], coeff)
            E_coeff = lambda a,b,c: scipy.integrate.dblquad(dE, 0,1, 0,1, args=(a,b,c), epsabs=0.0001, epsrel=0.0001)[0]
            compute = np.array(np.vectorize(E_coeff)(px,py,pz))
            E.append(compute.transpose())
        return E

    def plot_surface(self, ax):
        surf = sympy.lambdify([t, u], self.X.to_matrix(C))
        lin_t = np.linspace(0., 1., 10)
        lin_u = np.linspace(0., 1., 10)
        X, Y, Z = np.zeros(0), np.zeros(0), np.zeros(0)
        for dt in lin_t:
            for du in lin_u:
                o = surf(dt, du)
                X = np.append(X, o[0])
                Y = np.append(Y, o[1])
                Z = np.append(Z, o[2])
        ax.plot_trisurf(X,Y,Z)

system = System()
#point = Point(system, (0,0,5), 2*10**-6) # 2 microC
surface = Surface(system, t, u, 0)
#surface2 = Surface(system, t, u, 10)

#print("Point charge", surface.field_at_point(0,0,5))
#print(point.field_at_point(0,0,5))

#system.add_entity(point)
#system.add_entity(surface)
#system.add_entity(surface2)

point1 = Point(system, (0,4,0), 2E-6)
point2 = Point(system, (3,0,0), 10E-6)
system.add_entity(point1)
system.add_entity(point2)

force = point1.force_experienced()
print("Vector force", force)
print("Force =", np.sqrt(force.dot(force)))
