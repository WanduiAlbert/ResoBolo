<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

a = '0.00';
plottype = 'mag'

full_dir = '/twiki/pub/Resonators/20181008_NonLinearResonatorResponsivity/';

function update_fig_plots(){
    fig_link = full_dir+ "S21__"+plottype+"_plot_a"+a+".png";

    document["summary_fig"].src=fig_link;
}

function set_a(x){
    a = x;
    update_fig_plots();
}

function set_plottype(x){
    plottype = x;
    update_fig_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_a('0.00');">0.0</a></literal></sticky>|<sticky><literal><a href="javascript:set_a('0.50');">0.5</a></literal></sticky>|<sticky><literal><a href="javascript:set_a('0.90');">0.9</a></literal></sticky>|<sticky><literal><a href="javascript:set_a('1.50');">1.5</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_plottype('mag');">mag</a></literal></sticky>|<sticky><literal><a href="javascript:set_plottype('IQ');">IQ</a></literal></sticky>|<sticky><literal><a href="javascript:set_plottype('Re');">Re</a></literal></sticky>|<sticky><literal><a href="javascript:set_plottype('Im');">Im</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_fig;"> <img height=500  src="" name="summary_fig"></a> </literal></sticky>|

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_fig_plots();
-->
</script>
</literal></sticky>
