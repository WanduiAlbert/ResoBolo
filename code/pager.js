<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

result_freq = '298.6'

full_dir = '/twiki/pub/Resonators/20190910_NonEquilibriumQPGeneration/fig/';

function update_result_plots(){
    summary_result_flink = full_dir + "f_reso_fitted_" + result_freq + "MHz.png";

    document["summary_result_f"].src=summary_result_flink;

}

function set_result_freq(x){
    result_freq = x;
    update_result_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

| *Freq* |<sticky><literal><a href="javascript:set_result_freq('291.8');">291.8MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('298.6');">298.6MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('303.6');">303.6MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('317.9');">317.9MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('323.7');">323.7MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('332.5');">332.5MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('338.7');">338.7MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('354.8');">354.8MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('360.1');">360.1MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('367.6');">367.6MHz (R)</a></literal></sticky>|||||
| ^ |<sticky><literal><a href="javascript:set_result_freq('367.7');">367.7MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('372.3');">372.3MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('377.5');">377.5MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('388.7');">388.7MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('389.1');">389.1MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('394.8');">394.8MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('399.4');">399.4MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('412.7');">412.7MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('440.7');">440.7MHz (U)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('444.5');">444.5MHz (U)</a></literal></sticky>|||||
| ^ |<sticky><literal><a href="javascript:set_result_freq('449.9');">449.9MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('808.7');">808.7MHz (R)</a></literal></sticky>|<sticky><literal><a href="javascript:set_result_freq('905.5');">905.5MHz (U)</a></literal></sticky>|
|  |<sticky><literal> <a  href="javascript:location.href=summary_result_f;"> <img width=600  src="" name="summary_result_f"></a> </literal></sticky> ||||| <sticky><literal> <a  href="javascript:location.href=summary_result_Qi;"> <img width=600  src="" name="summary_result_Qi"></a> </literal></sticky> |||||||||

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_result_plots();
-->
</script>
</literal></sticky>

