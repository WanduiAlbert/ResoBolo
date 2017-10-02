# A set of class elements that represent the geometrical objects used in the 
# Caltech Intermediate Form for designing Integrated Circuits. The class definitions 
# used here will be used in writing scripts that make cif files.


class Point(object):
    """
    Represents a single coordinate location on a grid. Restricted to only using integers
    as the coordinate elements. If the coordinates are not specified while the Point is being 
    initialized, then the coordinates (0,0) are used.
    """
    def __init__(self, x=0, y=0):
        self.x = int(x)
        self.y = int(y)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)

    def __repr__(self):
        return "{0:d} {1:d}".format(self.x, self.y)

    def __str__(self):
        return "Point ({0:d}, {1:d})".format(self.x, self.y)


class Wire(object):
    """ 
    Represents a path on a mask.
    
    """
    def __init__(self, width, *points):
        """
        Wire(width, *points). Creates a wire object initialized to have a given width and passing 
        through the all the Point objects specified in the list points.
        """
        self.width = int(width)
        self.points = list(points)

    def addpoint(self, point):
        self.points.append(point)

    def end(self):
        return self.points[-1]

    def start(self):
        return self.points[0]

    def reverse(self):
        self.points = self.points[::-1]

    def __add__(self, other):
        if type(other) is type(self) and other.width != self.width:
            raise ValueError("Cannot combine two wires with different widths")
        total_points = self.points + other.points
        return Wire(self.width, *total_points)

    def __radd__(self, other):
        if type(other) is int:
            return self
        return self.__add__(other)

    def __repr__(self):
        start_str  = 'W {0:d} '.format(self.width)
        point_strs = map(repr, self.points)
        return (start_str + ' '.join(point_strs) + ';')

class Subroutine(object):
    """Represents a CIF subroutine that implements some particular functionality."""
    
    _start_code = 'DS'
    _end_code = 'DF'
    _cell_code = '9'

    def __init__(self, number, scale_num=0, scale_denom=0, cellname='', routine=''):
        self.number = int(number)
        self.scale_num = int(scale_num)
        self.scale_denom = int(scale_denom)
        self.cellname = cellname
        self.routine = routine
    
    def add_routine(self, routine):
        self.routine = routine

    def add_cellname(self, cellname):
        self.cellname=cellname

    def __repr__(self):
        ret_str = self._start_code + ' {0:d}'.format(self.number)
        if self.scale_num and  self.scale_denom:
            ret_str += ' {0:d} {1:d};\n'.format(self.scale_num, self.scale_denom)
        else:
            ret_str += ';\n'

        if self.cellname:
            ret_str += self._cell_code + ' ' + self.cellname + ';\n'

        # if not self.routine:
        #     raise SubRoutineMissingError('The implementation for this subroutine has not been given')

        ret_str += self.routine + '\n'
        ret_str += self._end_code + ';\n'

        return ret_str
        
class SubRoutineMissingError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
