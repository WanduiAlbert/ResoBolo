<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

fast_temperature = 'low'
fast_diagtemperature = 'low'
fast_channel = 'ch00'; 
fast_plot_type = 'nef';
fast_diagplot_type = 'meanpsd';
fast_txgain = '0'
 
full_dir = '/twiki/pub/Resonators/20180416_TLS_Screener_NoiseMeas/';

function update_figs_fastnoise(){
    summary_fastnoise_link = full_dir+ "fastchirp_fig/"+fast_temperature+"temp/fig/"+fast_plot_type+"_"+fast_channel+".png";

    document["summary_fastnoise"].src=summary_fastnoise_link;
}

function update_figs_diagfastnoise(){
    diag_fastnoise_link = full_dir+ "fastchirp_fig/"+fast_temperature+"temp/rawfig/noise_"+fast_txgain+"dB_"+fast_channel+"_"+fast_diagplot_type+".png";
    ts_fastnoise_link = full_dir+ "fastchirp_fig/"+fast_temperature+"temp/rawfig/noise_"+fast_txgain+"dB_ts.png";

    document["diag_fastnoise"].src=diag_fastnoise_link;
    document["ts_fastnoise"].src=ts_fastnoise_link;

}

function set_plottype_fastnoise(x){
    fast_plot_type = x;
    update_figs_fastnoise();
}

function set_channel_fastnoise(x){
    fast_channel = x;
    update_figs_fastnoise();
} 

function set_diagchannel_fastnoise(x){
    fast_channel = x;
    update_figs_diagfastnoise();
} 

function set_temperature_fastnoise(x){
    fast_temperature = x;
    update_figs_fastnoise();
}

function set_diagtemperature_fastnoise(x){
    fast_temperature = x;
    update_figs_diagfastnoise();
}

function set_txgain_fastnoise(x){
    fast_txgain = x;
    update_figs_diagfastnoise();
}

function set_diagplottype_fastnoise(x){
    fast_diagplot_type = x;
    update_figs_diagfastnoise();
}
// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_temperature_fastnoise('low');">100 mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_temperature_fastnoise('high');">2.8 K</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_plottype_fastnoise('nef');">nef</a></literal></sticky>|<sticky><literal><a href="javascript:set_plottype_fastnoise('nefasd');">nefasd</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_channel_fastnoise('ch00');">251.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_channel_fastnoise('ch01');">252.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_channel_fastnoise('ch02');">252.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_channel_fastnoise('ch03');">282.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_channel_fastnoise('ch04');">282.4MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_channel_fastnoise('ch05');">284.3MHz</a></literal></sticky>|


|  <sticky><literal> <a  href="javascript:location.href=summary_fastnoise;"> <img height=480  src="" name="summary_fastnoise"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_figs_fastnoise();
-->
</script>
</literal></sticky>

Diagnostic Plots

Noise Time Streams
</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_diagtemperature_fastnoise('low');">100 mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagtemperature_fastnoise('high');">2.8 K</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_txgain_fastnoise('0');">0 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('5');">5 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('10');">10 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('15');">15 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('20');">20 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('25');">25 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('30');">30 dB</a></literal></sticky>|

|  <sticky><literal> <a  href="javascript:location.href=ts_fastnoise;"> <img height=480  src="" name="ts_fastnoise"></a> </literal></sticky>  |

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_figs_diagfastnoise();
-->
</script>
</literal></sticky>

PSDs of the resonances
</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_diagtemperature_fastnoise('low');">100 mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagtemperature_fastnoise('high');">2.8 K</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_diagplottype_fastnoise('meanpsd');">mean asd</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagplottype_fastnoise('sampledpsd');">sampled psd</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_txgain_fastnoise('0');">0 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('5');">5 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('10');">10 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('15');">15 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('20');">20 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('25');">25 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_fastnoise('30');">30 dB</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_diagchannel_fastnoise('ch00');">251.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagchannel_fastnoise('ch01');">252.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagchannel_fastnoise('ch02');">252.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagchannel_fastnoise('ch03');">282.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagchannel_fastnoise('ch04');">282.4MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagchannel_fastnoise('ch05');">284.3MHz</a></literal></sticky>|


|  <sticky><literal> <a  href="javascript:location.href=diag_fastnoise;"> <img height=480  src="" name="diag_fastnoise"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_figs_diagfastnoise();
-->
</script>
</literal></sticky>


//Cross spectra matrix
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

cross_temperature = 'low'
cross_first_channel = 'ch01'; 
cross_second_channel = 'ch00'; 

// fast_plot_type = 'nef';
// fast_diagplot_type = 'meanpsd';
// fast_txgain = '0'
 
full_dir = '/twiki/pub/Resonators/20180416_TLS_Screener_NoiseMeas/';

function update_figs_cross_spectra(){
    summary_cross_link = full_dir+ "fastchirp_fig/"+cross_temperature+"temp/cross_spectra/rho_"+cross_first_channel+"_"+cross_second_channel+".png";

    document["summary_cross_spec"].src=summary_cross_link;
}

function set_cross_temperature(x){
    cross_temperature = x;
    update_figs_cross_spectra();
}

function set_first_channel(x){
    cross_first_channel = x;
    update_figs_cross_spectra();
} 

function set_second_channel(x){
    cross_second_channel = x;
    update_figs_cross_spectra();
} 

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_cross_temperature('low');">100 mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_cross_temperature('high');">2.8 K</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_first_channel('ch00');">251.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_first_channel('ch01');">252.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_first_channel('ch02');">252.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_first_channel('ch03');">282.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_first_channel('ch04');">282.4MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_first_channel('ch05');">284.3MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_second_channel('ch00');">251.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_second_channel('ch01');">252.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_second_channel('ch02');">252.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_second_channel('ch03');">282.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_second_channel('ch04');">282.4MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_second_channel('ch05');">284.3MHz</a></literal></sticky>|


|  <sticky><literal> <a  href="javascript:location.href=summary_cross_spec;"> <img height=480  src="" name="summary_cross_spec"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_figs_cross_spectra();
-->
</script>
</literal></sticky>

