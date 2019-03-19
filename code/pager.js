<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

timestream_voltage='1.388'
timestream_freq = '0.100'
 
full_dir = '/twiki/pub/Resonators/20190201_WaffleFab4TimeConstants/';

function update_timestream_plots(){
    summary_timestream_link = full_dir + timestream_voltage + "V_plots/rawfig/noise_"+timestream_voltage+"V_"+ timestream_freq+"Hz_ts.png";

    document["summary_timestream"].src=summary_timestream_link;

}

function set_timestream_voltage(x){
    timestream_voltage = x;
    update_timestream_plots();
}

function set_timestream_freq(x){
    timestream_freq = x;
    update_timestream_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

| *Bias Voltage* |<sticky><literal><a href="javascript:set_timestream_voltage('1.388');">1.388V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('1.963');">1.963V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('2.404');">2.404V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('2.776');">2.776V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('3.104');">3.104V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('3.400');">3.400V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('3.672');">3.672V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('3.926');">3.926V</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_voltage('4.164');">4.164V</a></literal></sticky>|

| *Freq* |<sticky><literal><a href="javascript:set_timestream_freq('0.100');">0.100Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.121');">0.121Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.146');">0.146Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.176');">0.176Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.212');">0.212Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.256');">0.256Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.309');">0.309Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.373');">0.373Hz</a></literal></sticky>|
|<sticky><literal><a href="javascript:set_timestream_freq('0.450');">0.450Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.543');">0.543Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.655');">0.655Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.791');">0.791Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('0.954');">0.954Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('1.151');">1.151Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('1.389');">1.389Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('1.677');">1.677Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('2.024');">2.024Hz</a></literal></sticky>|
|<sticky><literal><a href="javascript:set_timestream_freq('2.442');">2.442Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('2.947');">2.947Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('3.556');">3.556Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('4.292');">4.292Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('5.179');">5.179Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('6.251');">6.251Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('7.543');">7.543Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('9.103');">9.103Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('10.985');">10.985Hz</a></literal></sticky>|
|<sticky><literal><a href="javascript:set_timestream_freq('13.257');">13.257Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('15.999');">15.999Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('19.307');">19.307Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('23.300');">23.300Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('28.118');">28.118Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('33.932');">33.932Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('40.949');">40.949Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('49.417');">49.417Hz</a></literal></sticky>|<sticky><literal><a href="javascript:set_timestream_freq('59.636');">59.636Hz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_timestream;"> <img height=500  src="" name="summary_timestream"></a> </literal></sticky> |

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_timestream_plots();
-->
</script>
</literal></sticky>
