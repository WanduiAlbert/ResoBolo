import gdspy
import patches

"""
File: spreadgen_example.py
-----------------------------------------------------------------
An example of how to run the spreadsheet generation code in patches.py


"""







# File name to be read. Must be in GDS format.
fn = 'diplexer_FINAL_reformating.gds'
main_lib.read_gds(fn)
cell_dict = main_lib.cell_dict

# Give the name of the global overlay cell in the file just read.
globaloverlay = cell_dict['Module_features_just_tile']
# Specify all the mask files to be referenced as a list.
mask_list = [cell_dict['reticle_1'], cell_dict['reticle_2']]
# specify any cells in the global overlay that should be ignored because
# they are not on the reticles.
to_ignore = set(['frame'])
# A dictionary mapping the layer names to the GDSII layer numbers.
def_layers = {"Al_TES":1, "Outline":2, "PRO1":3, "Ti_TES":4, "GP":5, "PRO2":6,\
"VIA":7, "RES":8, "MS":9, "RES_pro":10, "heat_cap":11, "pic_frame":12,\
"LSN":13, "XeF2_STS":14, "New Layer":15, "Lens":20, "Die cutout":22,\
"Anodization":30}
# Specify the order in which the layers should be arranged in the
# spreadsheet.
layer_order = [1,3,4,6,5,8,10,7,9,13,11,12,14,30,2,15,20,22]
# This command creates a bunch of shots by running through the global
# overlay and the reticle list. If the cells on the mask are inverted
# specify so.
allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
		layer_dict=def_layers, layer_order=layer_order, cellsInverted=False)
# Create a patch table object from the list of shots generated.
patchtable = patches.PatchTable(allshots, 'diplexer_FINAL_reformating.xlsx')
patchtable.generate_spreadsheet()
