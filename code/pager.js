<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

freq = '300';

full_dir = '/twiki/pub/Resonators/20190109_SonnetCapacitorSimulations/';

function update_cap_plots(){
    summary_full_link = full_dir+ "Cap_"+freq+"MHz_with_boundaryY21_full.png";
    summary_fit_link = full_dir+ "Cap_"+freq+"MHz_with_boundaryY21.png";
    summary_res_link = full_dir+ "Cap_"+freq+"MHz_with_boundaryY21_residuals.png";

    document["summary_plot_full"].src=summary_full_link;
    document["summary_plot_fit"].src=summary_fit_link;
    document["summary_plot_res"].src=summary_res_link;
}

function set_freq(x){
    freq = x;
    update_cap_plots();
}


// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_freq('300');">300MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_freq('305');">305MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_freq('310');">310MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_freq('315');">315MHz</a></literal></sticky>| <sticky><literal><a href="javascript:set_freq('320');">320MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_freq('325');">325MHz</a></literal></sticky>| <sticky><literal><a href="javascript:set_freq('330');">330MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_freq('335');">335MHz</a></literal></sticky>| <sticky><literal><a href="javascript:set_freq('340');">340MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_freq('345');">345MHz</a></literal></sticky>|


|<sticky><literal> <a  href="javascript:location.href=summary_plot_full;"> <img height=500  src="" name="summary_plot_full"></a> </literal></sticky> |<sticky><literal> <a  href="javascript:location.href=summary_plot_fit;"> <img height=500  src="" name="summary_plot_fit"></a> </literal></sticky> | <sticky><literal> <a  href="javascript:location.href=summary_plot_res;"> <img height=500  src="" name="summary_plot_res"></a> </literal></sticky> |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_cap_plots();
-->
</script>
</literal></sticky>


