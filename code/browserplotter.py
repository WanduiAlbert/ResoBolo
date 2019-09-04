##############################
##                          ##
## RESONATOR VISUALIZER     ##
##                          ##
##############################

import numpy as np
from bokeh.layouts import gridplot, widgetbox
from bokeh.models import CustomJS, Slider
from bokeh.models.glyphs import Text
from bokeh.plotting import figure, ColumnDataSource
#from USRP_fitting import *
from bokeh.embed import components
from  yattag import Doc
from bokeh.models.widgets import RadioButtonGroup


class reso_plot(object):
    def __init__(self,resolution = 500, N_points = 1000, param_dict = None):
        self.N_points = N_points
        self.resolution = resolution
        if param_dict is None:
            self.param_dict = dict(
                f0 = 300,
                A=1,
                phi = 0,
                D = 0,
                Qi = 1e5,
                Qe = 1e5+0.j,
                Qr = 0,
                a = 0

            )
            self.param_dict['Qr'] = 1./(1./self.param_dict['Qi'] + 1./np.abs(self.param_dict["Qe"]))
        else:
            self.param_dict = param_dict

        # generate a frequency axis centered on f0
        freq = np.arange(self.N_points)*self.resolution + self.param_dict['f0']*1e6 - self.resolution*self.N_points/2
        # generate initial values for S21
        S21 = S21_func(freq, self.param_dict)

        real = S21.real
        imag = S21.imag
        mag_dB = 20*np.log10(np.abs(S21))
        phase = np.angle(S21)

        print("plotting...")

        self.source = ColumnDataSource(data=dict(
            real=real,
            imag=imag,
            mag_dB = mag_dB,
            phase = phase,
            freq = freq,
            x= [min(real) for i in range(len(freq))],
            y= [max(imag) for i in range(len(freq))],
            Qr = ["Qi = %.2fk"%(self.param_dict['Qr']/1e3)]+[None for i in range(len(freq)-1)]
        ))

        self.mag_plot = figure(plot_width=400, plot_height=400)#, output_backend="webgl"
        self.mag_plot.xaxis.axis_label = 'Frequency [Hz]'
        self.mag_plot.yaxis.axis_label = 'Magnitude [dB]'
        self.pha_plot = figure(plot_width=400, plot_height=400)
        self.pha_plot.xaxis.axis_label = 'Frequency [Hz]'
        self.pha_plot.yaxis.axis_label = 'Phase [rad]'
        self.iq_plot = figure(plot_width=400, plot_height=400, match_aspect=True)
        self.iq_plot.xaxis.axis_label = 'I [Arbitrary]'
        self.iq_plot.yaxis.axis_label = 'Q [Arbitrary]'

        self.mag_plot.line('freq', 'mag_dB', source=self.source, line_width=3, line_alpha=0.9)
        self.pha_plot.line('freq', 'phase', source=self.source, line_width=3, line_alpha=0.9)

        glyph = Text(x='x', y = 'y', text="Qr", angle=0, text_color="black")
        self.iq_plot.add_glyph(self.source, glyph)

        self.iq_plot.scatter('real','imag', source = self.source, size = 2)
        self.iq_plot.line('real', 'imag', source=self.source)

        callback = CustomJS(args=dict(source=self.source), code="""
            var Qi = Qi.value;
            var f0 = f0.value;
            var Qe_real = Qe_real.value;
            var Qe_imag = Qe_imag.value;
            var A = A.value;
            var a = a.value;
            var D = D.value;
            var phi = phi.value;
            
            var direction = dir.active;
            
            var Qr = 1.0/(1.0/Qi + 1.0/math.abs(math.complex(Qe_real,Qe_imag)))
            
            var data = source.data
            
            var f = data['freq'];
            var mag = data['mag_dB'];
            var phase = data['phase'];
            var I = data['imag'];
            var Q = data['real'];
            var Qe = math.complex(Qe_real, Qe_imag)
            var dQe = math.divide(1.0, Qe)
            
            complex_S21 = S21_func(f, f0, A, phi, D, math.divide(1,Qr), dQe.re, dQe.im, a, direction)
            
            for(var i = 0; i < f.length; i++){
                res_ =  math.subset(complex_S21, math.index(i))
                mag[i] = 20.*math.log10(math.abs(res_))
                phase[i] = math.atan2(res_.re,res_.im)
                //console.log("Val I: "+res_.re.toString()+" Q: "+res_.im.toString())
                I[i] = res_.re
                Q[i] = res_.im
            }
            
            
            var Qr_pos_y = data['y']
            var Qr_pos_x = data['x']
            var min_I = Math.max.apply(null, I)
            var min_Q = Math.min.apply(null, Q)
            for(var i = 0; i < f.length; i++){
                Qr_pos_y[i] = min_I
                Qr_pos_x[i] = min_Q
            }
            
            
            var Qr_text = data['Qr']
            Qr_text[0] =  "Qr: "+(Qr/1e3).toFixed(2)+"k"
            
            
            source.change.emit();
            
        """)

        f0_slider = Slider(start=max(min(freq)/1e6,1), end=max(freq)/1e6,
                           value=self.param_dict['f0'], step=1e-4, title="f0", callback=callback)
        callback.args["f0"] = f0_slider

        Qi_slider = Slider(start=max(self.param_dict['Qi']-1e5,1), end=self.param_dict['Qi']+1e5,
                           value=self.param_dict['Qi'], step=1e3,title="Qi", callback=callback)
        callback.args["Qi"] = Qi_slider

        Qe_real_slider = Slider(start=max(self.param_dict['Qe'].real - 1e5,1), end=self.param_dict['Qi']+1e5,
                           value=self.param_dict['Qe'].real, step=1e3, title="Qe real", callback=callback)
        callback.args["Qe_real"] = Qe_real_slider

        Qe_imag_slider = Slider(start=0, end=self.param_dict['Qi']+1e5,
                                value=self.param_dict['Qe'].imag, step=1e3, title="Qe imag", callback=callback)
        callback.args["Qe_imag"] = Qe_imag_slider

        A_slider = Slider(start=0.1, end=1,
                           value=self.param_dict['A'], step=1e-3, title="A", callback=callback)
        callback.args["A"] = A_slider

        phi_slider = Slider(start=0, end=np.pi,
                          value=self.param_dict['phi'], step=1e-3, title="phi", callback=callback)
        callback.args["phi"] = phi_slider

        D_slider = Slider(start=0, end=1,
                            value=self.param_dict['D'], step=1e-5, title="D", callback=callback)
        callback.args["D"] = D_slider

        a_slider = Slider(start=0, end=10,
                          value=self.param_dict['a'], step=1e-2, title="a", callback=callback)
        callback.args["a"] = a_slider

        radio_button_group = RadioButtonGroup(labels=['Low to high', 'High to low'], active=0, callback = callback)
        callback.args["dir"] = radio_button_group

        self.vbox = widgetbox(radio_button_group, f0_slider, Qi_slider, Qe_real_slider, Qe_imag_slider, A_slider, phi_slider, D_slider, a_slider)

        self.plots = {'mag': self.mag_plot, 'phase': self.pha_plot, 'iq': self.iq_plot, 'vbox': self.vbox}

        self.layout = gridplot(
            [self.mag_plot,self.pha_plot],
            [self.iq_plot,self.vbox]
        )




    def plot(self):



        script, div = components(self.layout)


        doc, tag, text = Doc().tagtext()

        doc.asis('<!DOCTYPE html>')
        with tag('html', lang="en"):

            with tag('head'):
                doc.asis('<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bokeh/0.13.0/bokeh.css" type="text/css" />')
                doc.asis(
                    '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bokeh/0.13.0/bokeh-widgets.min.css" />')
                with tag('script', src="https://cdnjs.cloudflare.com/ajax/libs/bokeh/0.13.0/bokeh.min.js"):
                    pass
                with tag('script', src="https://cdnjs.cloudflare.com/ajax/libs/bokeh/0.13.0/bokeh-widgets.js"):
                    pass
                with tag('script', src="https://unpkg.com/mathjs@5.4.2/dist/math.min.js"):
                    pass
                with tag('script', src="https://cdnjs.cloudflare.com/ajax/libs/bokeh/0.13.0/bokeh-gl.min.js"):
                    pass

                js_init =  """
                    var x0 = math.zeros(1)
                    var y0 = math.zeros(1)
                    var k1 = math.zeros(1)
                    var k2 = math.zeros(1)
                    var y_a =math.zeros(1)
                    var y_b =math.zeros(1)
                    var y1 = math.zeros(1)
                    var y2 = math.zeros(1)
                    var y = math.zeros(1)
                    //var y3 = math.zeros(1)
                    var x = math.zeros(1)
                    var S21_val = math.zeros(1)
                    var freq = []
                    var cable_phase = math.zeros(1)
                    
                    function S21_func(f, f0, A, phi, D, dQr, dQe_re, dQe_im, a, dir){
                        if (math.subset(math.size(x0),math.index(0)) != f.length){
                          console.log("Allocating variables: x0 is "+math.subset(math.size(x0),math.index(0)).toString()+" and f is: "+f.length.toString())
                          x0 = math.zeros(f.length)
                          y0 = math.zeros(f.length)
                          k1 = math.zeros(f.length)
                          k2 = math.zeros(f.length)
                          y_a = math.zeros(f.length)
                          y_b = math.zeros(f.length)
                          y1 = math.zeros(f.length)
                          y2 = math.zeros(f.length)
                          y = math.zeros(f.length)
                          //y3 = math.zeros(f.length)
                          x = math.zeros(f.length)
                          S21 = math.zeros(f.length)
                          cable_phase = math.zeros(f.length)
                          
                          for(var j = 0; j<f.length; j++){
                          freq.push(parseFloat(f[j]))
                          
                          }
                          //freq = math.zeros(f.length)
                          
                          console.log("Allocated.")

                          //console.log(freq)
                        }
                        
                        f0 = f0 * 1e6
                        dQe = math.complex(dQe_re , dQe_im)
                        
                        var cable_phase = math.exp(
                          math.multiply(
                            math.complex(0,2*Math.PI),
                            math.add(
                              math.multiply(1e-6 * D , math.subtract(freq, parseFloat(f0))),
                              phi
                            )
                          )
                        )
                    
                        x0 = math.dotDivide(math.subtract(freq , f0) , f0)
   
                        
                        y0 = math.dotDivide(x0 , dQr)
                        
                        k2 =  math.sqrt(math.subtract(
                                math.dotPow(
                                    math.add(
                                        math.dotDivide(
                                            math.dotPow(y0, 3.0 ),
                                            27.0
                                        ),
                                        math.dotDivide(y0, 12),
                                        parseFloat(a)/8.0
                                    ), 
                                    2.0
                                ),
                                math.dotPow(
                                    math.add(
                                        math.dotDivide(math.dotPow(y0, 2.0), 9.0),
                                        -1.0/12.0
                                    ),
                                    3.0
                                )
                        ))
                        
                        k1 = math.dotPow( math.add(parseFloat(a)/8.0 , math.dotDivide(y0,12.0), k2, math.dotDivide(math.dotPow(y0,3.0),27.0)), 1.0/3.0)
                        
                        
                        var eps = math.complex(-0.5,math.sqrt(3)/2.0)
                    
                        y_a = math.add(math.dotDivide(math.dotPow(y0,2.0),9.0),-1.0/12.0)
                        y_b = math.dotDivide(y0,3.0)
                        
                        y1 = math.add(y_b, math.dotDivide(y_a,k1), k1, math.complex(0,0))
                    
                        y2 = math.add(y_b, math.dotDivide(math.dotDivide(y_a,eps),k1), math.dotMultiply(eps,k1))
                    
                        //y3 = math.add(y_b, math.dotDivide(math.dotDivide(y_a,math.dotPow(eps,2.0)),k1), math.dotMultiply(math.dotPow(eps,2.0),k1) )
                        
                        
                        
                        var i;
                        for(i=0; i<f.length; i++){
                          //console.log(math.abs(k1[i]))
                          if(math.abs(math.subset(k1,math.index(i))) == 0.0){
                            y1["_data"][i] = y0["_data"][i]
                            y2["_data"][i] = y0["_data"][i]
                            //math.subset(y1, math.index(i), math.dotDivide(math.subset(y0,math.index(i)),3.0))
                            //math.subset(y2, math.index(i), math.dotDivide(math.subset(y0,math.index(i)),3.0))
                            //math.subset(y3, math.index(i), math.dotDivide(math.subset(y0,math.index(i)),3.0))
                          }
                        }
                    
                        
                        var low_to_high = dir;
                        
                        /* // got there with a radio button
                        if(f[0]<f[f.length-1]){
                          low_to_high = true
                        }else{
                          low_to_high = false
                        }
                        */
                        var thresh = 1e-4
                        
                        //var y = y1 //just initialization
                        
                        if (low_to_high == 0){
                          for(i=0; i<f.length; i++){
                            if (math.abs(math.subset(y2,math.index(i)).im) >= thresh){

                              y["_data"][i] = y1[i].re
                            }else{
                              y["_data"][i] = y2[i].re
                            }
                          }
                        }else{
                          for(i=0; i<f.length; i++){
                            if (math.abs(math.subset(y1,math.index(i)).im) >= thresh){
                              y["_data"][i] = y2[i].re
                            }else{
                              y["_data"][i] = y1[i].re
                            }
                          }
                        }
                        
                        x = math.dotMultiply(y, dQr)
                        
                        
                        
                        S21_val =  math.dotMultiply(
                            parseFloat(A),
                            math.dotMultiply( 
                                cable_phase,
                                math.subtract(
                                    1.0,
                                    math.dotDivide(
                                        dQe,
                                        math.add(
                                            dQr,
                                            math.dotMultiply(math.complex(0,2.0), x)
                                        )
                                    )
                                )
                            )
                        )

                        return S21_val
                    }
                """
                with tag('script'):
                    doc.asis(js_init)
                doc.asis(script)

            with tag('body'):

                doc.asis(div)



        with open("result_page.html", "w") as f:
            f.write(doc.getvalue())


if __name__ == "__main__":
    x = reso_plot()
    x.plot()
