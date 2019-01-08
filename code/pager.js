<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

tc_freq = '428.7';
tc_power = '2.5';

full_dir = '/twiki/pub/Resonators/20181105_TKIDWaffleTimeConstantMeas/';

function update_timeconstant_plots(){
    summary_tc_link = full_dir+ tc_freq+"MHz_plots/fig/reso_"+tc_freq+"MHz_noise_25dB_power"+tc_power+"pW_timeconstant.png";

    document["summary_tc_plot"].src=summary_tc_link;
}

function set_tc_freq(x){
    tc_freq = x;
    update_timeconstant_plots();
}

function set_tc_power(x){
    tc_power = x;
    update_timeconstant_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_tc_freq('428.7');">428.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_tc_freq('469.9');">469.9MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_tc_power('2.5');">2.5pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_tc_power('5.0');">5.0pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_tc_power('7.5');">7.5pW</a></literal></sticky>|<sticky><literal><a href="javascript:set_tc_power('10.0');">10.0pW</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_tc_plot;"> <img height=500  src="" name="summary_tc_plot"></a> </literal></sticky> |

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_timeconstant_plots();
-->
</script>
</literal></sticky>


<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

frequency = '428.7';
 
vfull_dir = '/twiki/pub/Resonators/20181002_WaffleTKIDGMeasurements/';

function update_lowrawfig_plots(){
    summary_lowrawfig_link = vfull_dir+ "lowtemp/rawfig/reso_"+frequency+"MHz.png";

    document["summary_lowrawfig"].src=summary_lowrawfig_link;
}

function set_frequency(x){
    frequency = x;
    update_lowrawfig_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_frequency('428.7');">428.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('469.9');">469.9MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_lowrawfig;"> <img height=800  src="" name="summary_lowrawfig"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_lowrawfig_plots();
-->
</script>
</literal></sticky>
