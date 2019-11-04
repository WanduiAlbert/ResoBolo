
---+++ VNA Scans

<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

frequency = '275.40';
 
vfull_dir = '/twiki/pub/Resonators/20191023_CF69CNetanal/';

function update_vna_plots(){
    summary_vna_link = vfull_dir+ "powersweep/reso_"+frequency+"MHz.png";

    document["summary_vna"].src=summary_vna_link;
}

function set_frequency(x){
    frequency = x;
    update_vna_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_frequency('275.40');">275.40MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('279.80');"> 279.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('296.70');">296.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('302.90');">302.90MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('314.20');">314.20MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('324.00');">324.00MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('689.80');">689.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('954.50');">954.50MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_frequency('1194.30');">1194.30MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1311.80');">1311.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1316.30');">1316.30MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1352.50');">1352.50MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1457.60');">1457.60MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1465.20');">1465.20MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1568.80');">1568.80MHz</a></literal>|<sticky><literal><a href="javascript:set_frequency('1614.00');">1614.00MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_frequency('1646.90');">1646.90MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1718.70');">1718.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1757.60');">1757.60MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('1897.80');">1897.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('450.30');">450.30MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_vna;"> <img height=600  src="" name="summary_vna"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_plots();
-->
</script>
</literal></sticky>

---+++ diagnostic fit plots
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

diag_frequency = '275.40';
diag_txgain = '-110'
 
dfull_dir = '/twiki/pub/Resonators/20191023_CF69CNetanal/';

function update_vna_diagnostic_plots(){
    summary_vna_diag_link = dfull_dir+ "powersweep/diagnostic/reso_"+diag_frequency+"MHz_"+diag_txgain+"dBm.png";

    document["summary_vna_diag"].src=summary_vna_diag_link;
}

function set_diag_frequency(x){
    diag_frequency = x;
    update_vna_diagnostic_plots();
}

function set_diag_txgain(x){
    diag_txgain = x;
    update_vna_diagnostic_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>
|<sticky><literal><a href="javascript:set_diag_frequency('275.40');">275.40MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('279.80');"> 279.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('296.70');">296.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('302.90');">302.90MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('314.20');">314.20MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('324.00');">324.00MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('689.80');">689.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('954.50');">954.50MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_diag_frequency('1194.30');">1194.30MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1311.80');">1311.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1316.30');">1316.30MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1352.50');">1352.50MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1457.60');">1457.60MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1465.20');">1465.20MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1568.80');">1568.80MHz</a></literal>|<sticky><literal><a href="javascript:set_diag_frequency('1614.00');">1614.00MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_diag_frequency('1646.90');">1646.90MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1718.70');">1718.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1757.60');">1757.60MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('1897.80');">1897.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('450.30');">450.30MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_diag_txgain('-110');">-110dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-108');">-108dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-106');">-106dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-104');">-104dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-102');">-102dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-100');">-100dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-98');">-98dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-96');">-96dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-94');">94dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-92');">-92dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-90');">-90dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-88');">-88dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-86');">-86dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-84');">-84dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-82');">-82dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_txgain('-80');">-80dBm</a></literal></sticky> | 

|<sticky><literal> <a  href="javascript:location.href=summary_vna_diag;"> <img height=600  src="" name="summary_vna_diag"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_diagnostic_plots();
-->
</script>
</literal></sticky>

---++ Temperature sweeps
warm attenuation: 10dB cold attenuation: 60dB Power out = -20 dBm

data: ~/data/20190603_darkresosweep/tempsweep/autonetanal_20190603_092821.txt

---+++ Summary Plots
   * Tc_hist.png: <br />
     <img src="%ATTACHURLPATH%/Tc_hist.png" alt="Tc_hist.png" height="500" />

   * frvsTc.png: <br />
     <img src="%ATTACHURLPATH%/frvsTc.png" alt="frvsTc.png" height="500" />

   * all_resonators_QivsT_fitted.png: <br />
     <img src="%ATTACHURLPATH%/all_resonators_QivsT_fitted.png" alt="all_resonators_QivsT_fitted.png" height="500" />

   * reso_124_to_514mK_dfvsT.png: <br />
     <img src="%ATTACHURLPATH%/reso_124_to_514mK_dfvsT.png" alt="reso_124_to_514mK_dfvsT.png" height="500" />

   * reso_124_to_514mK_QivsT.png: <br />
     <img src="%ATTACHURLPATH%/reso_124_to_514mK_QivsT.png" alt="reso_124_to_514mK_QivsT.png" height="500" />

---+++ Temperature sweep VNAs

<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

temp_frequency = '275.40';
 
vfull_dir = '/twiki/pub/Resonators/20191023_CF69CNetanal/';

function update_temp_vna_plots(){
    summary_temp_vna_link = vfull_dir+ "tempsweep/rawfig/reso_"+temp_frequency+"MHz.png";

    document["summary_temp_vna"].src=summary_temp_vna_link;
}

function set_temp_frequency(x){
    temp_frequency = x;
    update_temp_vna_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_temp_frequency('275.40');">275.40MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('279.80');"> 279.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('296.70');">296.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('302.90');">302.90MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('314.20');">314.20MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('324.00');">324.00MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('689.80');">689.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('954.50');">954.50MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_temp_frequency('1194.30');">1194.30MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1311.80');">1311.80MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1316.30');">1316.30MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1352.50');">1352.50MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1457.60');">1457.60MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1465.20');">1465.20MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1568.80');">1568.80MHz</a></literal>|<sticky><literal><a href="javascript:set_temp_frequency('1614.00');">1614.00MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_temp_frequency('1646.90');">1646.90MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1718.70');">1718.70MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1757.60');">1757.60MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('1897.80');">1897.80MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_temp_vna;"> <img height=600  src="" name="summary_temp_vna"></a> </literal></sticky>  |

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_temp_vna_plots();
-->
</script>
</literal></sticky>

---+++ Fits to the Temp Sweep

I limited the data used in this analysis to T < 360 mK since the fits get spurious for very broad resonances. I fit for both alphak and Tc and use MCMC to marginalize over the alpha distribution to get best estimate for Tc. The marginalized distributions are pretty tightly constrained even with all the degeneracy between alphak and Tc.

<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

fit_x_ch = '0';

fit_x_full_dir = '/twiki/pub/Resonators/20191023_CF69CNetanal/';

function update_fit_x_plots(){
    summary_fit_x_link = fit_x_full_dir+ "tempsweep/fig/reso_"+fit_x_ch+"_xvsT_fitted.png";
    summary_fit_qi_link = fit_x_full_dir+ "tempsweep/fig/reso_"+fit_x_ch+"_QivsT_fitted.png";

    document["summary_fit_x"].src=summary_fit_x_link;
    document["summary_fit_qi"].src=summary_fit_qi_link;

}

function set_fit_x_ch(x){
    fit_x_ch = x;
    update_fit_x_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_fit_x_ch('0');">275.4MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('1');"> 279.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('2');">296.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('3');">302.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('4');">314.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('5');">324.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('6');">689.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('7');">954.5MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_fit_x_ch('8');">1194.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('9');">1311.8MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('10');">1316.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('11');">1352.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('12');">1457.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('13');">1465.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('14');">1568.8MHz</a></literal>|<sticky><literal><a href="javascript:set_fit_x_ch('15');">1614.0MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_fit_x_ch('16');">1646.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('17');">1718.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('18');">1757.6MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_fit_x_ch('19');">1897.8MHz</a></literal></sticky>|

| <sticky><literal> <a  href="javascript:location.href=summary_fit_x;"> <img height=600  src="" name="summary_fit_x"></a> </literal></sticky> | <sticky><literal> <a  href="javascript:location.href=summary_fit_qi;"> <img height=600  src="" name="summary_fit_qi"></a> </literal></sticky>  |



<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_fit_x_plots();
-->
</script>
</literal></sticky>
