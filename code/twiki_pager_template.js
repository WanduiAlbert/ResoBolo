<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

frequency_resonator = '253.3'
power_resonator = '-110'
power_range_resonator = 'low'

power_name_rawspectra = '0.00pW';
txgain_name_rawspectra = '15dB';
ch_name_rawspectra = 'all';
measurement_name_rawspectra = 'gain';

full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/';

function update_figs_netanals(){
    summary_netanals_link = full_dir+power_range_resonator+"fig/"+"_reso_"+frequency_resonator+"MHz.png";
    diag_netanals_link = full_dir+"fig/"+power_range_resonator+"_reso_"+frequency_resonator+"MHz_"+power_resonator+"dBm.png";
    document["summary_netanals"].src=summary_netanals_link;
    document["diag_netanals"].src=diag_netanals_link;
}

function set_range_netanals(x){
    power_range_resonator = x;
    update_figs_netanals();
}

function set_frequency_netanals(x){
    frequency_resonator = x;
    update_figs_netanals();
}

function set_power_netanals(x){
    power_resonator = x;
    update_figs_netanals();
}


function update_figs_rawspectra(){
    rawspectra_link = full_dir+"fig/"+measurement_name_rawspectra+"_"+txgain_name_rawspectra+"_"+power_name_rawspectra+"_"+ch_name_rawspectra+"_meanpsd.png";
    document["rawspectra"].src=rawspectra_link;
}

// function set_power_rawspectra(x) {
//     power_name_rawspectra = x;
//     update_figs_rawspectra();
// }
// function set_txgain_rawspectra(x) {
//     txgain_name_rawspectra = x;
//     update_figs_rawspectra();
// }
// function set_ch_rawspectra(x) {
//     ch_name_rawspectra = x;
//     update_figs_rawspectra();
// }
// function set_measurement_rawspectra(x) {
//     measurement_name_rawspectra = x;
//     update_figs_rawspectra();
// }

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_range_netanals('low');">low</a></literal></sticky>|<sticky><literal><a href="javascript:set_range_netanals('high');">high</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_frequency_netanals('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_netanals('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_netanals('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_netanals('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_netanals('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_netanals('285.9');">285.9MHz</a></literal></sticky>|

// |<sticky><literal><a href="javascript:set('all');">all</a></literal></sticky>|<sticky><literal><a href="javascript:set_ch_rawspectra('ch00');">ch00</a></literal></sticky>|<sticky><literal><a href="javascript:set_ch_rawspectra('ch01');">ch01</a></literal></sticky>|<sticky><literal><a href="javascript:set_ch_rawspectra('ch02');">ch02</a></literal></sticky>|<sticky><literal><a href="javascript:set_ch_rawspectra('ch03');">ch03</a></literal></sticky>|<sticky><literal><a href="javascript:set_ch_rawspectra('ch04');">ch04</a></literal></sticky>|

// |<sticky><literal><a href="javascript:set_measurement_rawspectra('gain');">gain</a></literal></sticky>|<sticky><literal><a href="javascript:set_measurement_rawspectra('noise');">noise</a></literal></sticky>|

|  <sticky><literal> <a  href="javascript:location.href=summary_netanals;"> <img height=480  src="" name="summary_netanals"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_figs_rawspectra();
-->
</script>
</literal></sticky>



// Template for the noise data
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

frequency = '251.6';
plot_type = 'nef';
diag_plot_type = 'waterfall';
txgain = '0';
 
 

full_dir = '/twiki/pub/Resonators/20180416_TLS_Screener_NoiseMeas/';

function update_figs_noise(){
    summary_noise_link = full_dir+ "fig/reso_"+frequency+"MHz_"+plot_type+".png";
    diag_noise_link = full_dir+"fig/diagnostics/"+diag_plot_type+"_"+frequency+"MHz_0Vdc_"+txgain+"dB.png";
    document["summary_noise"].src=summary_noise_link;
    document["diag_noise"].src=diag_noise_link;
}

function set_plottype_noise(x){
    plot_type = x;
    update_figs_noise();
}

function set_frequency_noise(x){
    frequency = x;
    update_figs_noise();
}

function set_txgain_noise(x){
    txgain = x;
    update_figs_noise();
}

function set_diagplottype_noise(x){
    diag_plot_type = x;
    update_figs_noise();
} 

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_plottype_noise('nef');">nef</a></literal></sticky>|<sticky><literal><a href="javascript:set_plottype_noise('nefasdonfloor');">nefasdonfloor</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_frequency_noise('251.6');">251.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('252.3');">252.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('253.1');">253.1MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('281.7');">281.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('282.6');">282.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('284.3');">284.3MHz</a></literal></sticky>|


|  <sticky><literal> <a  href="javascript:location.href=summary_noise;"> <img height=480  src="" name="summary_noise"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_figs_noise();
-->
</script>
</literal></sticky>


</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_diagplottype_noise('waterfall');">waterfall</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagplottype_noise('section');">section</a></literal></sticky>|<sticky><literal><a href="javascript:set_diagplottype_noise('ts');">timestream</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_frequency_noise('251.6');">251.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('252.3');">252.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('253.1');">253.1MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('281.7');">281.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('282.6');">282.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency_noise('284.3');">284.3MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_txgain_noise('0');">0 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_noise('5');">5 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_noise('10');">10 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_noise('15');">15 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_noise('20');">20 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_noise('25');">25 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain_noise('30');">30dB</a></literal></sticky>|


|  <sticky><literal> <a  href="javascript:location.href=diag_noise;"> <img height=480  src="" name="diag_noise"></a> </literal></sticky>  |


<sticky></td></tr></table>