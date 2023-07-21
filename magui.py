import tkinter

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np
import sympy
from sympy import sin, cos, pi, tan, log

import Electromagnet

from importlib import reload
Electromagnet = reload(Electromagnet)

def make_system():
    system = Electromagnet.System()

    point1 = Electromagnet.Point(system, (0,4,0), 2E-6)
    point2 = Electromagnet.Point(system, (3,0,0), 10E-6)

    #system.add_entity(point1)
    #system.add_entity(point2)

    t,u = sympy.symbols("t u")

    plane1 = Electromagnet.Surface(system, t-0.5, u-0.5, -1)
    plane2 = Electromagnet.Surface(system, t-0.5, u-0.5, 1, charge_density=-1)

    system.add_entity(plane1)
    system.add_entity(plane2)

    return system

def line_system():
    t = sympy.symbols("t")
    system = Electromagnet.System()
    #line = Electromagnet.Line(system, 0,0,t, t_limits=(-5, 5), current=1, charge_density=1)
    line = Electromagnet.Line(system, 0,0.8*t,10*t**3, t_limits=(-0.5, 0.5), current=1, charge_density=1)
    system.add_entity(line)
    return system

root = tkinter.Tk()
root.wm_title("Embedding in Tk")

fig = Figure()
ax = fig.add_subplot(projection="3d")
t = np.arange(0, 3, .01)
#line, = ax.plot(t, 2 * np.sin(2 * np.pi * t))
ax.set_xlabel("x / m")
ax.set_ylabel("y / m")
ax.set_zlabel("z / m")

E_system = make_system()
#E_system = line_system()
E_system.plot_entities(ax)
E_system.plot_electric_field(ax)
#E_system.plot_magnetic_field(ax)

ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

#ax.plot(0.8,1, marker="o")

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()

# pack_toolbar=False will make it easier to use a layout manager later on.
toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
toolbar.update()

canvas.mpl_connect(
    "key_press_event", lambda event: print(f"you pressed {event.key}"))
canvas.mpl_connect("key_press_event", key_press_handler)


# Packing order is important. Widgets are processed sequentially and if there
# is no space left, because the window is too small, they are not displayed.
# The canvas is rather flexible in its size, so we pack it last which makes
# sure the UI controls are displayed as long as possible.
toolbar.pack(side=tkinter.BOTTOM, fill=tkinter.X)
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

tkinter.mainloop()
