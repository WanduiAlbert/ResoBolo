<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

rho_txgain = '0';
rho_power = '0.0';
rho_channel1 = '00';
rho_channel2 = '02';

full_dir = '/twiki/pub/Resonators/20181009_WaffleTKIDNEPMeasurements/';

function update_rho_plots(){
	rho_fig_link = full_dir+"cross_spectra/cm_coeff_noise_"+rho_txgain+"dB_"+rho_power+"pW_ch"+rho_channel1+"_ch"+rho_channel2+".png";
	csd_fig_link = full_dir+"cross_spectra/cm_noise_"+rho_txgain+"dB_"+rho_power+"pW_ch"+rho_channel1+"_ch"+rho_channel2+".png";

    document["summary_csd_main"].src=csd_fig_link;
    document["summary_rho_main"].src=rho_fig_link;
}

function set_rho_txgain(x){
    rho_txgain = x;
    update_rho_plots();
}

function set_rho_power(x){
    rho_power = x;
    update_rho_plots();
}

function set_rho_channel1(x){
    rho_channel1 = x;
    update_rho_plots();
}

function set_rho_channel2(x){
    rho_channel2 = x;
    update_rho_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>


|<sticky><literal><a href="javascript:set_rho_txgain('0');">0dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_txgain('5');">5dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_txgain('10');">10dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_txgain('15');">15dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_txgain('20');">20dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_txgain('25');">25dB</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_rho_power('0.0');">0.0pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_power('2.5');">2.5pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_power('5.0');">5.0pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_power('7.5');">7.5pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_power('10.0');">10.0pW</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_rho_channel1('00');">428.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_channel1('01');">449.10MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_channel1('02');">469.90MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_rho_channel2('00');">428.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_channel2('01');">449.10MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_rho_channel2('02');">469.90MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_rho_main;"> <img height=500  src="" name="summary_rho_main"></a> </literal></sticky>| <sticky><literal> <a  href="javascript:location.href=summary_csd_main;"> <img height=500  src="" name="summary_csd_main"></a> </literal></sticky>|

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_rho_plots();
-->
</script>
</literal></sticky>
