<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

summary_index = '1';

vfull_dir = '/twiki/pub/Resonators/20190104_DarkResoArrayMeasTc/';

function update_summary_vna_plots(){
    summary_vna_link_x = vfull_dir+ "fig/reso_"+summary_index+"_xvsT_fitted.png";
    summary_vna_link_Q = vfull_dir+ "fig/reso_"+summary_index+"_QivsT_fitted.png";

    document["summary_vna_x"].src=summary_vna_link_x;
    document["summary_vna_Q"].src=summary_vna_link_Q;
}

function set_summary_index(x){
    summary_index = x;
    update_summary_vna_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_summary_index('1');">290.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('2');">303.1MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('3');">310.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('4');">312.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('5');">328.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('6');">363.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('7');">367.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('8');">382.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('9');">386.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('10');">391.1MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('11');">400.4MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('12');">448.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_summary_index('13');">458.5MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_summary_vna_x;"> <img height=600  src="" name="summary_summary_vna_x"></a> </literal></sticky>  | <sticky><literal> <a  href="javascript:location.href=summary_summary_vna_Q;"> <img height=600  src="" name="summary_summary_vna_Q"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_summary_vna_plots();
-->
</script>
</literal></sticky>

