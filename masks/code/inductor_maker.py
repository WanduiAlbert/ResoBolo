#! /usr/bin/env python3
from cif import Point, Wire, Subroutine
import operator
from functools import reduce

def make_inductor(start, ind_length, ind_width, wire_width, wire_spacing):
    inductor = Wire(wire_width, start)

    # My strategy will be to create as many long wires as I can and then join them at the
    # ends using shorter wire segments to complete the inductor.

    # Use integer division to obtain the actual number of wires. Integer division gives
    # the number of spacings. Add 1 to get the number of long connectors
    num_long = (ind_width // (wire_width + wire_spacing)) + 1
    N_wires  = 2*num_long - 1 

    wires = [0]* N_wires
    curr_start = start
    for i in range(0, N_wires, 2):
        wires[i] = Wire(wire_width, curr_start, curr_start + Point(ind_length, 0))
        curr_start += Point(0, wire_spacing + wire_width)

    # Now we need to join the wires at the end in an alternating fashion

    # Join left wires
    for i in range(3, N_wires, 4):
        wires[i] = Wire(wire_width, wires[i-1].start(), wires[i+1].start())

    # Join right wires
    for i in range(1, N_wires, 4):
        wires[i] = Wire(wire_width, wires[i-1].end(), wires[i+1].end())
        wires[i + 1].reverse()

    return sum(wires)


if __name__ == "__main__":
    start_coord_lower = Point(15771000, 3010500)
    ind_length = 301000
    ind_width_lower = 55000
    wire_width = 1000
    wire_spacing = 1000
    inductor_lower = make_inductor(start_coord_lower, ind_length, ind_width_lower, wire_width, wire_spacing)

    start_coord_upper = inductor_lower.end() + Point(0, 21000)
    ind_width_upper = 54000
    inductor_upper = make_inductor(start_coord_upper, ind_length, ind_width_upper, wire_width, wire_spacing)
    
    # We will join the two inductors with a single wire connection
    ind_connector = Wire(wire_width, inductor_lower.end(), inductor_upper.start())

    inductor = inductor_lower  + ind_connector + inductor_upper

    base_routine = Subroutine(1, 1, 10, 'INDUCTOR', 'L ALUMINUM;\n' +  repr(inductor))


    with open('inductor_coded.cif', 'w') as f:
        f.write('(Inductor coded using a python script)\n')
        f.write(repr(base_routine))
        f.write('E;')