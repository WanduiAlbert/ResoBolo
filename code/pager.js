<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

fastchirpres_ch = '00';
fastchirpres_tx = '0';
fastchirpres_power = '0.00';

fcfull_dir = '/twiki/pub/Resonators/20181009_WaffleTKIDNEPMeasurements/';

function update_fastchirpresnoise_plots(){
    summary_fastchirpres_link = fcfull_dir+ "rawfig/noise_"+fastchirpres_tx+"dB_"+fastchirpres_power+"pW_ch"+fastchirpres_ch+"_meanpsd_fitted.png";

    document["summary_fastchirpres_plot"].src=summary_fastchirpres_link;
}

function set_fastchirpres_channel(x){
    fastchirpres_ch = x;
    update_fastchirpresnoise_plots();
}

function set_fastchirpres_tx(x){
    fastchirpres_tx = x;
    update_fastchirpresnoise_plots();
}

function set_fastchirpres_power(x){
    fastchirpres_power = x;
    update_fastchirpresnoise_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_fastchirpres_channel('00');">428.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_channel('01');">449.1MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_channel('02');">469.9MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_fastchirpres_tx('0');">0dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_tx('5');">5dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_tx('10');">10dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_tx('15');">15dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_tx('20');">20dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_tx('25');">25dB</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_fastchirpres_power('0.00');">0.0pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_power('2.50');">2.5pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_power('5.00');">5.0pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_power('7.50');">7.5pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_fastchirpres_power('10.00');">10.0pW</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_fastchirpres_plot;"> <img height=600  src="" name="summary_fastchirpres_plot"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_fastchirpresnoise_plots();
-->
</script>
</literal></sticky>
