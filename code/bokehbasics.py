from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure, output_file, show

x = [1, 2, 3, 4, 5]
y = [6, 8, 9, 16, 1]

output_file("lines.html", title='Dumb example')

p = figure(title="Simple Line Plot", x_axis_label='x', y_axis_label='y')

p.line(x, y, legend_label='Temp.', line_width=2)

show(p)
