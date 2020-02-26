

<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

maxtemp = '100.0';

vfull_dir = '/twiki/pub/Resonators/20191108_AlphaTcDegeneracy/';

function update_max_plots(){
    summary_chisq_link = vfull_dir+ "chisq/chisq_vs_alphaTc_Tmax"+maxtemp+"mK.png";
    summary_dQchisq_link = vfull_dir+ "chisq/chisq_vs_alphaTc_Tmax"+maxtemp+"mK.png";

    document["summary_chisq"].src=summary_chisq_link;
    document["summary_dQchisq"].src=summary_dQchisq_link;
}

function set_maxtemp(x){
    maxtemp = x;
    update_chisq_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_maxtemp('100.0');">100.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('150.0');"> 150.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('200.0');">200.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('250.0');">250.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('300.0');">300.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('350.0');">350.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('400.0');">400.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('450.0');">450.0mK</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_maxtemp('500.0');">500.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('550.0');">550.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('600.0');">600.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('650.0');">650.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('700.0');">700.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('750.0');">750.0mK</a></literal></sticky>|<sticky><literal><a href="javascript:set_maxtemp('800.0');">800.0mK</a></literal>|

|<sticky><literal> <a  href="javascript:location.href=summary_chisq;"> <img height=600  src="" name="summary_chisq"></a> </literal></sticky>  | <sticky><literal> <a  href="javascript:location.href=summary_dQchisq;"> <img height=600  src="" name="summary_dQchisq"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_chisq_plots();
-->
</script>
</literal></sticky>

