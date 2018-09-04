import gdspy
import numpy as np
import pdb
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, Border, Side, PatternFill
from string import ascii_uppercase

allowed_element_types = set([gdspy.Cell, gdspy.CellReference, gdspy.CellArray])
default = np.zeros(2)
scale = 1000

class cellNode():

    """
    Keeps track of a single cell and its origin relative to the root
    of the cell hierarchy tree
    """
    def __init__(self, cell, origin=default):
        self.cell = cell
        self.origin = origin

    def update_origin(self, new_origin):
        self.origin += new_origin

# scales all the values in an array by the scale
def scalearr(arr, scale):
  return np.asarray(arr)/scale


# Need this functionality to work with excel spreadsheets
chars = list(ascii_uppercase)
chars2 = [x + y for x in chars for y in chars]
chars += chars2
char2num = {char:(ind+1) for ind, char in enumerate(chars)}
num2char = {(ind+1):char for ind, char in enumerate(chars)}

al = Alignment(horizontal="center", wrap_text=False)
font = Font(name='Times New Roman', size=10)
thin = Side(border_style='thin', color='000000')
border = Border(top=thin, left=thin, right=thin, bottom=thin)
fill = PatternFill(fill_type=None, start_color='FFFFFFFF', end_color='FF000000')

#def as_text(value): return str(value) if value is not None else ""

def style_range(ws, cell_range, border=Border(), fill=None, font=None, alignment=None):
    """
    :param ws:  Excel worksheet instance
    :param range: An excel range to style (e.g. A1:F20)
    :param border: An openpyxl Border
    :param fill: An openpyxl PatternFill or GradientFill
    :param font: An openpyxl Font object
    """



    rows = ws[cell_range]
    dims = {}
    for row in rows:
        for cell in row:
            cell.font = font
            cell.border = border
            cell.fill = fill
            cell.alignment = alignment
    #         if cell.value:
    #             dims[cell.column] = max((dims.get(cell.column, 0), len(cell.value)))
    # for col, value in dims.items():
    #     ws.column_dimensions[col].width = value


    #first_cell = ws[cell_range.split(":")[0]]
    #if alignment:
    #    ws.merge_cells(cell_range)
    #    first_cell.alignment = alignment

    #rows = ws[cell_range]
    #if font:
    #    first_cell.font = font

    #for cell in rows[0]:
    #    cell.border = cell.border + top
    #for cell in rows[-1]:
    #    cell.border = cell.border + bottom

    #for row in rows:
    #    l = row[0]
    #    r = row[-1]
    #    l.border = l.border + left
    #    r.border = r.border + right
    #    if fill:
    #        for c in row:
    #            c.fill = fill

def populate_column(sheet, col, startrow, dataset):
    N = len(dataset)
    for index in range(N):
        row = startrow + index
        sheet.cell(column = col, row=row, value = dataset[index])


class PatchTable():

    # Some class constants that define the relative positions of the TVPA and
    # AGA Marks with respect to the center of the alignment marks cell
    # These are defined for alignment_marks_patch_new_r
    #tvpa_shift = np.array([0.125, 0])
    #aga_x_shift = np.array([-0.125, 0])
    #aga_y_shift = np.array([-0.225, 0])

    # For alignment_marks_patch_r
    tvpa_shift = np.array([0.16875, 0])
    aga_x_shift = np.array([-0.08125, 0])
    aga_y_shift = np.array([-0.18125, 0])


    def __init__(self, shotlist, wb_filename="ResonatorArray.xlsx"):
        self.shots = shotlist
        self.filename = wb_filename
        self.wb = Workbook()
        self.layers = np.array(list(map(lambda x: x.get_layer(),\
                self.shots)))
        self.names = np.array(list(map(lambda x: x.get_name(),\
                self.shots)))
        self.mask_shifts = np.array(list(map(lambda x: x.get_mask_shift(),\
                self.shots)))
        self.cell_sizes = np.array(list(map(lambda x: x.get_cell_size(),\
                self.shots)))
        self.cell_shifts = np.array(list(map(lambda x: x.get_cell_shift(),\
                self.shots)))
        self.cell_bboxs = np.array(list(map(lambda x: x.get_cell_bbox(),\
                self.shots)))
        self.mask_names = np.array(list(map(lambda x: x.get_maskname(),\
                self.shots)))
        self.array_sizes = np.array(list(map(lambda x: x.get_array_size(),\
                self.shots)))
        self.array_centers = np.array(list(map(lambda x:\
                x.get_array_center(), self.shots)))
        self.array_stepsizes = np.array(list(map(lambda x:\
                x.get_array_stepsize(), self.shots)))
        self.array_shifted = np.array(list(map(lambda x:\
                x.get_array_shifts(), self.shots)))
        self.array_positions = np.array(list(map(lambda x:\
                x.get_array_positions(), self.shots)))
        self.calculated_xys = np.array(list(map(lambda x:\
                x.get_calculated_xy(), self.shots)))
        self.desired_xys = np.array(list(map(lambda x:\
                x.get_desired_xy(), self.shots)))
        self.shifts_xys = np.array(list(map(lambda x:\
                x.get_shift_xy(), self.shots)))
        self.num_shots = len(self.layers)
        self.xl = self.mask_shifts[:, 0] + self.cell_bboxs[:,0,0]
        self.xr = self.mask_shifts[:, 0] + self.cell_bboxs[:,1,0]
        self.yu = self.mask_shifts[:, 1] + self.cell_bboxs[:,1,1]
        self.yd = self.mask_shifts[:, 1] + self.cell_bboxs[:,0,1]
        self.align_layer = self.layers[0]

        self.wafer_shifts = self.mask_shifts - self.cell_shifts
        self.ws = self.wb.active
        self.ws.title = 'Patches' #get the active worksheet in which to save the data
        self.shift_ws = self.wb.create_sheet('Shift_Calculator')

    def populate_alignment_info(self):
        self.ws['B4'] = 'AGA Marks'
        self.ws['H4'] = 'TVPA'
        self.ws['B5'] = 'direction'
        self.ws['C5'] = 'layer'
        self.ws['D5'] = 'x'
        self.ws['E5'] = 'y'
        self.ws['G5'] = 'layer'
        self.ws['H5'] = 'x'
        self.ws['I5'] = 'y'
        self.ws['B6'] = 'y'
        self.ws['B7'] = 'x'
        self.ws['C6'] = self.align_layer
        self.ws['C7'] = self.align_layer
        self.ws['G6'] = self.align_layer

        # Get the relative position of the AGA and TVPA marks from the cell
        # shift of the alignment layer cell
        align_shift = self.cell_shifts[0]
        tvpa_pos = align_shift + PatchTable.tvpa_shift
        aga_x_pos = align_shift + PatchTable.aga_x_shift
        aga_y_pos = align_shift + PatchTable.aga_y_shift

        self.ws['D6'] = aga_y_pos[0]
        self.ws['E6'] = aga_y_pos[1]
        self.ws['D7'] = aga_x_pos[0]
        self.ws['E7'] = aga_x_pos[1]
        self.ws['H6'] = tvpa_pos[0]
        self.ws['I6'] = tvpa_pos[1]

    def populate_shiftcalc(self):
        ws = self.shift_ws

        cell_range = "A1:AA200"
        style_range(ws, cell_range, border=border, fill=fill,\
                font=font, alignment=al)
        ws['A2'] = 'Job'
        ws['B2'] = 'patch'
        ws['C1'] = 'Array'
        ws['E1'] = 'Centered (mm)'
        ws['G1'] = 'Pitch (mm)'
        ws['J1'] = 'Column'
        ws['J2'] = 'position'
        ws['O1'] = 'Row'
        ws['O2'] = 'position'
        ws['K1'] = 'Calculate'
        ws['P1'] = 'Calculate'
        ws['L1'] = 'Desired'
        ws['Q1'] = 'Desired'
        ws['M1'] = 'Shift in'
        ws['R1'] = 'Shift in'
        ws['C2'] = 'C'
        ws['D2'] = 'R'
        ws['E2'] = 'x'
        ws['G2'] = 'x'
        ws['K2'] = 'x'
        ws['L2'] = 'x'
        ws['M2'] = 'x'
        ws['F2'] = 'y'
        ws['H2'] = 'y'
        ws['P2'] = 'y'
        ws['Q2'] = 'y'
        ws['R2'] = 'y'
        N = len(self.array_shifted)
        start_row = 3
        row = start_row
        for i in range(N):
            if not self.array_shifted[i]:
                continue
            ws['A' + str(row)] = self.names[i]
            ws['B' + str(row)] = self.layers[i]
            ws['C' + str(row)] = self.array_sizes[:,0][i]
            ws['D' + str(row)] = self.array_sizes[:,1][i]
            ws['E' + str(row)] = self.array_centers[:, 0][i]
            ws['F' + str(row)] = self.array_centers[:, 1][i]
            ws['G' + str(row)] = self.array_stepsizes[:, 0][i]
            ws['H' + str(row)] = self.array_stepsizes[:, 1][i]

            Ncols = len(self.array_positions[:, 0][i])
            Nrows = len(self.array_positions[:, 1][i])
            #print (self.calculated_xys[:, 0][i])
            #print (self.calculated_xys[:, 1][i])
            crow = row
            for j in range(Ncols):
                ws['J' + str(crow)] = self.array_positions[:, 0][i][j]
                ws['K' + str(crow)] = self.calculated_xys[:, 0][i][j]
                ws['L' + str(crow)] = self.desired_xys[:, 0][i][j]
                ws['M' + str(crow)] = self.shifts_xys[:, 0][i][j]
                crow += 1
            rrow = row
            for j in range(Nrows):
                ws['O' + str(rrow)] = self.array_positions[:, 1][i][j]
                ws['P' + str(rrow)] = self.calculated_xys[:, 1][i][j]
                ws['Q' + str(rrow)] = self.desired_xys[:, 1][i][j]
                ws['R' + str(rrow)] = self.shifts_xys[:, 1][i][j]
                rrow += 1
            row += max([Ncols, Nrows]) + 2
            #row += 1
        ws.merge_cells('C1:D1')
        ws.merge_cells('E1:F1')
        ws.merge_cells('G1:H1')
        ws.merge_cells('E1:F1')

        # Adjust the spreadsheet to ensure that the width of the cells is enough
        for column_cells in ws.columns:
            length = max(len(str(cell.value or ""))+2 for cell in column_cells)
            if length < 6: length = 6
            ws.column_dimensions[column_cells[0].column].width = length



    def generate_spreadsheet(self):
        # Write down the TVPA and AGA information at the top of the spreadsheet
        self.populate_alignment_info()

        # Format the spreadsheet before you write the data. Formatting multiple
        # cells seems to overwrite the data
        self.startrow = 15
        self.endrow = 15 + self.num_shots

        cell_range = "A4:AB%d"%(self.endrow)
        style_range(self.ws, cell_range, border=border, fill=fill,\
                font=font, alignment=al)

        populate_column(self.ws, char2num['A'], 15, self.layers)
        populate_column(self.ws, char2num['C'], 15, self.names)
        populate_column(self.ws, char2num['D'], 15, self.mask_shifts[:,0])
        populate_column(self.ws, char2num['E'], 15, self.mask_shifts[:,1])
        populate_column(self.ws, char2num['F'], 15, self.cell_sizes[:, 0])
        populate_column(self.ws, char2num['G'], 15, self.cell_sizes[:, 1])
        populate_column(self.ws, char2num['M'], 15, self.cell_shifts[:, 0])
        populate_column(self.ws, char2num['N'], 15, self.cell_shifts[:, 1])

        populate_column(self.ws, char2num['H'], 15, self.xl)
        populate_column(self.ws, char2num['I'], 15, self.xr)
        populate_column(self.ws, char2num['J'], 15, self.yu)
        populate_column(self.ws, char2num['K'], 15, self.yd)

        populate_column(self.ws, char2num['O'], 15, self.wafer_shifts[:,0])
        populate_column(self.ws, char2num['P'], 15, self.wafer_shifts[:,1])
        populate_column(self.ws, char2num['S'], 15, self.mask_names)

        populate_column(self.ws, char2num['T'], 15, self.array_sizes[:,0])
        populate_column(self.ws, char2num['U'], 15, self.array_sizes[:,1])
        populate_column(self.ws, char2num['V'], 15, self.array_centers[:, 0])
        populate_column(self.ws, char2num['W'], 15, self.array_centers[:, 1])
        populate_column(self.ws, char2num['X'], 15, self.array_stepsizes[:, 0])
        populate_column(self.ws, char2num['Y'], 15, self.array_stepsizes[:, 1])
        populate_column(self.ws, char2num['AB'], 15, self.array_shifted)

        # Adjust the spreadsheet to ensure that the width of the cells is enough
        for column_cells in self.ws.columns:
            length = max(len(str(cell.value or ""))+2 for cell in column_cells)
            if length < 6: length = 6
            self.ws.column_dimensions[column_cells[0].column].width = length
        
        self.ws['A14'] = 'Layer'
        self.ws['B14'] = 'Description'
        self.ws['C14'] = 'Cell name'
        self.ws['D13'] = 'Mask Shift'
        self.ws['D14'] = 'x'
        self.ws['E14'] = 'y'
        self.ws['F13'] = 'Cell size'
        self.ws['F14'] = 'x'
        self.ws['G14'] = 'y'
        self.ws['M13'] = 'Cell Shift'
        self.ws['M14'] = 'x'
        self.ws['N14'] = 'y'
        self.ws['H13'] = 'Blade Coordinates'
        self.ws['H14'] = 'xl'
        self.ws['I14'] = 'xr'
        self.ws['J14'] = 'yu'
        self.ws['K14'] = 'yd'
        self.ws['O13'] = 'Wafer Shift'
        self.ws['O14'] = 'x'
        self.ws['P14'] = 'y'
        self.ws['Q14'] = 'Patch'
        self.ws['Q15'] = 'Frontside'
        self.ws['R14'] = 'Dose'
        self.ws['S14'] = 'Mask'
        self.ws['T12'] = 'Alternate Array Layout'
        self.ws['T13'] = 'size'
        self.ws['T14'] = 'C'
        self.ws['U14'] = 'R'
        self.ws['V13'] = 'center'
        self.ws['V14'] = 'x'
        self.ws['W14'] = 'y'
        self.ws['X13'] = 'step size'
        self.ws['X14'] = 'x'
        self.ws['Y14'] = 'y'
        self.ws['Z12'] = 'To Skip'
        self.ws['Z13'] = '(1-indexed, from LL)'
        self.ws['Z14'] = 'C'
        self.ws['AA14'] = 'R'

        self.ws.merge_cells('D13:E13')
        self.ws.merge_cells('F13:G13')
        self.ws.merge_cells('M13:N13')
        self.ws.merge_cells('H13:K13')
        self.ws.merge_cells('O13:P13')
        self.ws.merge_cells('T12:Y12')
        self.ws.merge_cells('T13:U13')
        self.ws.merge_cells('V13:W13')
        self.ws.merge_cells('X13:Y13')
        self.ws.merge_cells('Z12:AA12')
        self.ws.merge_cells('Z13:AA13')

        self.populate_shiftcalc()

        self.wb.save(self.filename)

class Shot():
    def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,\
            "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9,\
            "XeF2":10, "GP":12, 'Wafer Outline':22}
    inv_layers = {value:key for key, value in def_layers.items()}
    layer_list = sorted(list(def_layers.values()))
    index = list(range(len(layer_list)))
    ordering = {x:y for x, y in zip(layer_list, index)}

    def update_layers(layers_dict):
        Shot.def_layers = layers_dict
        Shot.inv_layers = {value:key for key, value in Shot.def_layers.items()}
        Shot.layer_list = sorted(list(Shot.def_layers.values()))
        Shot.index = list(range(len(Shot.layer_list)))
        Shot.ordering = {x:y for x, y in zip(Shot.layer_list, Shot.index)}

    def update_layerorder(layer_order):
        Shot.index = list(range(len(layer_order)))
        Shot.ordering = {x:y for x, y in zip(layer_order, Shot.index)}


    def __init__(self, cell, cell_shift, cell_size, cell_bbox, mask_name="", isArray=False, **kwargs):
        try:
            assert(type(cell) == gdspy.Cell or type(cell) == gdspy.CellReference)
        except AssertionError:
            print ("Please provide a valid cell to construct the shot")
        self.cell = cell
        self.cellname = cell.name
        layers = list(cell.get_layers())
        if len(layers) > 1 :
          raise RuntimeError ("Cannot construct a shot from a cell with multiple layers")
        self.layer = layers[0]
        self.cell_bbox = np.asarray(cell_bbox)
        self.cell_size = np.asarray(cell_size)
        self.mask_name = mask_name
        self.cell_shift = np.asarray(cell_shift)
        self.isArray = isArray
        #self.blade_coords = {'xl':, 'xr':, 'yt':, 'yb':}
        if isArray:
            self.ncols = kwargs['num_cols']
            self.nrows = kwargs['num_rows']
            self.center = np.asarray(kwargs['center'])
            self.xspacing = kwargs['xspacing']
            self.yspacing = kwargs['yspacing']
            self.is_shifted = False
            try:
                if kwargs['is_shifted']:
                    self.is_shifted = True
                    self.column_pos = kwargs['column_pos']
                    self.row_pos = kwargs['row_pos']
                    self.calculated_x = kwargs['calculated_x']
                    self.calculated_y = kwargs['calculated_y']
                    self.desired_x = kwargs['desired_x']
                    self.desired_y = kwargs['desired_y']
                    self.shift_x = kwargs['shift_x']
                    self.shift_y = kwargs['shift_y']
            except KeyError:
                #do nothing
                print ('Nope. Nothing here')
        # defaults
        self.maskcell=""
        self.maskcellsize=default
        self.mask_shift=default

    def update_mask_location(self, maskcellref, maskname):
        self.mask_name = maskname
        maskcellname = maskcellref.ref_cell.name
        assert(maskcellname.startswith(self.cellname))

        self.maskcell = maskcellref.ref_cell
        self.maskcellsize = scalearr(get_size(maskcellref), scale)
        self.mask_shift = scalearr(maskcellref.origin, scale)

    def __repr__(self):
        if not self.isArray:
            return "Shot({0}, {1}, {2}, {3})".format(self.cellname,\
                    self.cell_shift,self.mask_shift, self.maskcellsize)
        return "Shot({0}, {1}, {2}, {3}, {4}, {5})".format(self.cellname,\
            self.center, self.mask_shift, self.maskcellsize,\
            (self.ncols, self.nrows),(self.xspacing, self.yspacing))

    def __eq__(self, other):
        if Shot.ordering[self.layer] == Shot.ordering[other.layer]:
            return self.cellname == other.cellname
        return False

    def __lt__(self, other):
        if Shot.ordering[self.layer] == Shot.ordering[other.layer]:
            if self.cellname.startswith('alignment'): return True
            if other.cellname.startswith('alignment'): return False
            return self.cellname < other.cellname
        return Shot.ordering[self.layer] < Shot.ordering[other.layer]

    def get_layer(self):
        return Shot.inv_layers[self.layer]

    def get_name(self):
        return self.cellname

    def get_maskname(self):
        return self.mask_name

    def get_mask_shift(self):
        return self.mask_shift

    def get_cell_size(self):
        return self.maskcellsize

    def get_cell_shift(self):
        return self.cell_shift

    def get_cell_bbox(self):
        return self.cell_bbox

    def get_array_size(self):
        if not self.isArray: return None, None
        return self.ncols, self.nrows

    def get_array_center(self):
        if not self.isArray: return None, None
        return self.center

    def get_array_stepsize(self):
        if not self.isArray: return None, None
        return self.xspacing, self.yspacing

    def get_array_shifts(self):
        if not self.is_shifted:
            return ""
        else:
            return "shift"

    def get_array_positions(self):
        if not self.is_shifted: return None, None
        return self.column_pos, self.row_pos

    def get_calculated_xy(self):
        if not self.is_shifted: return None, None
        return self.calculated_x, self.calculated_y

    def get_desired_xy(self):
        if not self.is_shifted: return None, None
        return self.desired_x, self.desired_y

    def get_shift_xy(self):
        if not self.is_shifted: return None, None
        return self.shift_x, self.shift_y

def translate(shot, shift):
    shot.cell_shift -= np.asarray(shift)

def get_size(cell):
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    dx = (xmax - xmin)
    dy = (ymax - ymin)
    return np.array([dx, dy])

def get_center(cell):
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    x0 = (xmax + xmin)//2
    y0 = (ymax + ymin)//2
    return np.array([x0, y0])

def same_mask_cellname(shot):
    return shot.cellname

def inv_mask_cellname(shot):
    if shot.cellname.endswith('_r'):
        name = shot.cellname
    elif shot.cellname.endswith('inv'):
        name = shot.cellname
    elif shot.cellname.endswith('-r'):
        name = shot.cellname
    else:
        name = shot.cellname + '_r'
    return name

def gen_patches_table(globaloverlay, mask_list, ignored_cells, layer_dict=None,\
        layer_order=None, cellsInverted=True):
    Shot.update_layers(layer_dict)
    if layer_order:
        Shot.update_layerorder(layer_order)
    gcomponents = globaloverlay.elements
    allshots = []
    # pdb.set_trace()
    for component in gcomponents:
        if type(component) not in allowed_element_types: continue
        if component.ref_cell.name in ignored_cells: continue
        shots_made = makeshot(component, mask_list=mask_list)
        #print (component.ref_cell.name, len(shots_made))
        allshots.extend(shots_made)

    mcomponents = {}
    for mask in mask_list:
        mcomponents[mask.name] = [x for x in mask.elements if type(x) in allowed_element_types]

    for shot in allshots:
        if cellsInverted:
            name = inv_mask_cellname(shot)
        else:
            name = same_mask_cellname(shot)
        for mask in mcomponents:
            try:
                match = list(filter(lambda x: x.ref_cell.name == name,
                    mcomponents[mask]))[0]
                if not match: continue
                shot.update_mask_location(match, mask)
            except IndexError:
                # print ("Could not find a matching cell on mask for cell {:s}".format(name))
                continue

    allshots.sort()
    return allshots

empty_dict = dict()

def get_cell_asymmetry(cell):
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    dx = xmax - np.abs(xmin)
    dy = ymax - np.abs(ymin)
    return [dx, dy]

def get_array_shifts(element, parent_args):
    # Unpack the parent array info for easier computation
    Nc,Nr = parent_args['num_cols'], parent_args['num_rows']
    X0, Y0 = parent_args['center']
    DX, DY = parent_args['xspacing'], parent_args['yspacing']

    # Get info of the child array
    nc, nr = element.columns, element.rows
    x0, y0 = scalearr(get_center(element) - get_center(element.ref_cell), scale)
    dx, dy = scalearr(element.spacing, scale)

    # Compute the new array properties
    Mc, Mr = Nc*nc, Nr*nr
    Mx = np.arange(Mc) + 1
    My = np.arange(Mr) + 1
    nx = (Mx - 1)%nc + 1
    Nx = (Mx - nx)/nc + 1
    ny = (My - 1)%nr + 1
    Ny = (My - ny)/nr + 1

    calcx = X0 - dx/2*(Mc - 2*Mx + 1)
    calcy = Y0 + dy/2*(Mr - 2*My + 1)
    shiftx = -(Nc - 2*Nx + 1)*(DX - nc*dx)/2
    shifty = -(Nr - 2*Ny + 1)*(DY - nr*dy)/2
    desiredx = calcx + shiftx
    desiredy = calcy + shifty

    shiftArray = True
    # I need to be careful of the cases where either dx or dy is zero
    if nc == 1: dx = 0
    if nr == 1: dy = 0

    if dx == 0:
        calcx = np.zeros_like(Mx)
        shiftx = np.zeros_like(Mx)
        desiredx = np.zeros_like(Mx)
        dx = DX

    if dy == 0:
        calcy = np.zeros_like(My)
        shifty = np.zeros_like(My)
        desiredy = np.zeros_like(My)
        dy = DY

    if np.all(shiftx) == 0:
        calcx = np.zeros_like(Mx)
        shiftx = np.zeros_like(Mx)
        desiredx = np.zeros_like(Mx)
        DX = dx

    if np.all(shifty) == 0:
        calcy = np.zeros_like(My)
        shifty = np.zeros_like(My)
        desiredy = np.zeros_like(My)
        DY = dy

    # Sometimes nested arrays end up falling on a regular grid with exactly no
    # row/column shifts. In this case, we will intuit the existence of this new array and
    # use it.
    if dx == DX and dy == DY:
        shiftArray = False

    new_args = {'center':parent_args['center'], 'num_cols':Mc, 'num_rows':Mr,\
            'xspacing':dx, 'yspacing':dy, 'is_shifted':shiftArray,\
        'column_pos':Mx, 'calculated_x':calcx, 'desired_x':desiredx,\
        'shift_x': shiftx, 'row_pos':My, 'calculated_y':calcy,\
        'desired_y':desiredy, 'shift_y': shifty}
    return new_args

def makeshot(curr_element, parent_origin=default, parentIsArray=False, arrayArgs=empty_dict, mask_list=[]):
    #if curr_element.ref_cell.name == "DetBiasLeft": pdb.set_trace()
    #if curr_element.ref_cell.name in ignored_cells: return
    if type(curr_element) not in allowed_element_types:
        return []
    curr_cell = curr_element.ref_cell
    cell_center = get_center(curr_cell)
    curr_origin = curr_element.origin
    abs_origin = parent_origin + scalearr(curr_origin, scale)
    cell_shift = abs_origin
    cell_size = scalearr(get_size(curr_cell), scale)
    cell_bbox = scalearr(curr_cell.get_bounding_box(), scale)

    isArray = False
    if type(curr_element) == gdspy.CellArray and not parentIsArray:
        # Need to correct the array center by the diff btn the dimensions of the ref_cell
        arr_center = get_center(curr_element)
        arr_center = scalearr(arr_center - cell_center , scale) + parent_origin
        abs_origin = default#scalearr(cell_center, scale) # Really crucial. All subcells of the cell shot
        cell_shift = abs_origin
        # in the array have their cell shifts relative to this cell.
        xspacing, yspacing = scalearr(curr_element.spacing, scale)
        args = {'num_cols':curr_element.columns, 'num_rows':curr_element.rows,\
            'center':arr_center, 'xspacing':xspacing, 'yspacing':yspacing,
            'is_shifted':False}
        isArray = True

    elif type(curr_element) == gdspy.CellReference and parentIsArray:
        args = arrayArgs
        isArray = True

    elif type(curr_element) == gdspy.CellReference and not parentIsArray:
        args = {'num_cols':1, 'num_rows':1, 'center':cell_shift, 'xspacing':0,\
            'yspacing':0, 'is_shifted':False}
        cell_shift = default

    elif type(curr_element) == gdspy.CellArray and parentIsArray:
        abs_origin = default
        arr_center = scalearr(get_center(curr_element), scale)
        #if curr_element.ref_cell.name == "L_MS_feed_bolometer":
        #    print (arr_center)
        cell_shift = arr_center
        arr_center = scalearr(arr_center - cell_center , scale) + parent_origin
        args = get_array_shifts(curr_element, arrayArgs)
        isArray = True


    haschildren = bool(curr_cell.get_dependencies())

    # Sometimes cells have children but are on the mask.  We want to
    # use the largest mask cell possible, so we stop digging early
    # if we find a cell on the mask
    cellonmask = False
    for m in mask_list:
        if curr_cell in m.get_dependencies():
            cellonmask = True

    if cellonmask or (not haschildren):

        try:
            # if curr_element.ref_cell.name == 'gnd_lf' or curr_element.ref_cell.name == 'nitride_lf':
            #     pdb.set_trace()
            shot = Shot(curr_cell, cell_shift, cell_size, cell_bbox, isArray=True, **args)
            return [shot]
        except RuntimeError:
            print ("Failed for", curr_element)
            return []

    child_elements = curr_cell.elements
    # If a cell has dependencies, then the elements gives a list of all the cell
    # references in this cell. If the cell has no dependencies, then subelements
    # just gives the polygon set that makes up the cell. I don't want the
    # polygonset stuff.
    child_shotlist = []
    for child in child_elements:
        if type(child) not in allowed_element_types: continue
        # If the current shot is part of a larger array, I want to keep track of
        # that information
        child_shotlist.extend(makeshot(child, abs_origin, isArray, args, mask_list=mask_list))
    return child_shotlist


# def makeshot(element, parent=None, hierarchy=0):
#     #pdb.set_trace()
#     isArray = False
#     if type(element) == gdspy.CellArray:
#         isArray = True
#     elif type(element) != gdspy.CellReference:
#         return []

#     cell = element.ref_cell
#     cellname = cell.name
#     # If a cell has dependencies, then the elements gives a list of all the cell
#     # references in this cell. If the cell has no dependencies, then subelements
#     # just gives the polygon set that makes up the cell. I don't want the
#     # polygonset stuff.
#     subelements = []
#     if cell.get_dependencies():
#         subelements = cell.elements

#     if isArray:
#         arr_center = scalearr(get_center(element), scale)
#         ncols = element.columns
#         nrows = element.rows
#         xspacing, yspacing = scalearr(element.spacing, scale)
#         args = {'num_cols':ncols, 'num_rows':nrows,\
#             'center':arr_center, 'xspacing':xspacing, 'yspacing':yspacing}

#     # If deps is an empty set then immediately construct the shot from the
#     # element. Otherwise loop through the dependencies of the current element and
#     # obtain shots from each of them.
#     if not subelements:
#         cell_shift = scalearr(element.origin, scale)
#         cell_size = scalearr(get_size(element), scale)
#     else:
#         shotlist = []
#         for el in subelements:
#             shots = makeshot(el, parent=element, hierarchy=hierarchy + 1)
#             # For each shot that has been made from each element, update its origin
#             # relative to the origin for  the element
#             if hierarchy >= 2:
#                 translate(shot, scalearr(parent.origin, scale))
#             shotlist.extend(shots)
#         # If the current shot is part of a larger array, I want to keep track of
#         # that information
#         if isArray:
#             for shot in shotlist:
#                 shot.isArray = True
#                 shot.center = arr_center
#                 shot.ncols = ncols
#                 shot.nrows = nrows
#                 shot.xspacing, shot.yspacing = xspacing, yspacing
#         return shotlist

#     if isArray:
#         args = {'num_cols':ncols, 'num_rows':nrows,\
#             'center':arr_center, 'xspacing':xspacing, 'yspacing':yspacing}
#         shot = Shot(cell, cell_shift, cell_size, isArray=True, **args)
#     else:
#         args = {'num_cols':1, 'num_rows':1,\
#             'center':cell_shift, 'xspacing':0, 'yspacing':0}
#         try:
#             shot = Shot(cell, (0, 0), cell_size, isArray=True, **args)
#         except RuntimeError:
#             print ("Failed for", element)
#             return []

#     if hierarchy >= 2:
#         translate(shot, scalearr(parent.origin, scale))
#     return [shot]















