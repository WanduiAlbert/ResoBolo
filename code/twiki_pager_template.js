<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

frequency = '320';
txgain = '0';

sfull_dir = '/twiki/pub/Resonators/20180605_TLS_Screener_SingleTone/vna_noise_comparison/';

function update_vna_noise_comp_plots(){
    summary_vna_link = full_dir+ "VNA_vs_Noise_"+frequency+"MHz_txgain_"+ txgain+ "dB.png";

    document["summary_vna"].src=summary_vna_link;
}

function set_frequency(x){
    frequency = x;
    update_vna_noise_comp_plots();
}

function set_txgain(x){
    txgain = x;
    update_vna_noise_comp_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_frequency('320');">320MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('360');">360MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_txgain('0');">0 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain('5');">5 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain('10');">10 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain('15');">15 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain('20');">20 dB</a></literal></sticky>|<sticky><literal><a href="javascript:set_txgain('25');">25 dB</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_vna;"> <img height=480  src="" name="summary_vna"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_qp_single_plots();
-->
</script>
</literal></sticky>


---++++ VNA scans
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

frequency = '253.3';
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/low_power/';

function update_vna_plots(){
    summary_vna_link = full_dir+ "rawfig/reso_"+frequency+"MHz.png";

    document["summary_vna"].src=summary_vna_link;
}

function set_frequency(x){
    frequency = x;
    update_vna_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_frequency('394.3');">394.3MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_vna;"> <img height=480  src="" name="summary_vna"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_plots();
-->
</script>
</literal></sticky>

---++++ diagnostic fit plots
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

diag_frequency = '253.3';
diag_txgain = '-110'
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/low_power/';

function update_vna_diagnostic_plots(){
    summary_vna_diag_link = full_dir+ "rawfig/diagnostic/reso_"+diag_frequency+"MHz_"+diag_txgain+"dBm.png";

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

|<sticky><literal><a href="javascript:set_diag_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_frequency('394.3');">394.3MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_diag_txgain('-110');">-110dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-108');">-108dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-106');">106dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-104');">104dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-102');">102dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-100');">100dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-98');">98dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-96');">96dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-94');">94dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-92');">92dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-90');">90dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-88');">88dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-86');">86dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_txgain('-84');">84dBm</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_vna_diag;"> <img height=480  src="" name="summary_vna_diag"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_diagnostic_plots();
-->
</script>
</literal></sticky>


---+++ -90 to -74 dBm measurements
Measurements with 20dB warm attenuation + 60dB cold attenuation: ~/data/20180413_netanal/autonetanal_20180413_125709.txt

TLS effects seen more clearly. 

| <img src="%ATTACHURLPATH%/reso_-90_to_-74dBm_QivsP.png" alt="reso_-90_to_-74dBm_QivsP.png" width="620" /> | <img src="%ATTACHURLPATH%/reso_-90_to_-74dBm_fvsP.png" alt="reso_-90_to_-74dBm_fvsP.png" width="550" /> | 

---++++ Individual Q vs P plots
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

single_high_frequency = '254.5';
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/high_power/';

function update_qp_high_single_plots(){
    summary_qp_high_single_link = full_dir+ "fig/diagnostics/resonator_"+single_high_frequency+"MHz_QvsP.png";

    document["summary_qp_high_single"].src=summary_qp_high_single_link;
}

function set_single_high_frequency(x){
    single_high_frequency = x;
    update_qp_high_single_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_single_high_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_single_high_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_single_high_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_single_high_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_single_high_frequency('284.0');">284.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_single_high_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_single_high_frequency('394.3');">394.3MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_qp_high_single;"> <img height=480  src="" name="summary_qp_high_single"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_qp_high_single_plots();
-->
</script>
</literal></sticky>


---++++ VNA scans
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

high_frequency = '253.3';
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/high_power/';

function update_vna_high_plots(){
    summary_vna_high_link = full_dir+ "rawfig/reso_"+high_frequency+"MHz.png";

    document["summary_vna_high"].src=summary_vna_high_link;
}

function set_high_frequency(x){
    high_frequency = x;
    update_vna_high_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_high_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_high_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_high_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_high_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_high_frequency('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_high_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_high_frequency('394.3');">394.3MHz</a></literal></sticky>|

|<sticky><literal> <a  href="javascript:location.href=summary_vna_high;"> <img height=480  src="" name="summary_vna_high"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_high_plots();
-->
</script>
</literal></sticky>

---++++ diagnostic fit plots
<literal>
<script type=''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''text/javascript''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''>// <![CDATA[

diag_high_frequency = '253.3';
diag_high_txgain = '-90'
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/high_power/';

function update_vna_high_diagnostic_plots(){
    summary_vna_high_diag_link = full_dir+ "rawfig/diagnostic/reso_"+diag_high_frequency+"MHz_"+diag_high_txgain+"dBm.png";

    document["summary_vna_high_diag"].src=summary_vna_high_diag_link;
}

function set_diag_high_frequency(x){
    diag_high_frequency = x;
    update_vna_high_diagnostic_plots();
}

function set_diag_high_txgain(x){
    diag_high_txgain = x;
    update_vna_high_diagnostic_plots();
}

// ]]></script>

</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_diag_high_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_frequency('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_frequency('394.3');">394.3MHz</a></literal></sticky>|

|<sticky><literal><a href="javascript:set_diag_high_txgain('-90');">90dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_txgain('-88');">88dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_txgain('-86');">86dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_txgain('-84');">84dBm</a></literal></sticky>|<sticky><literal><a href="javascript:set_diag_high_txgain('-82');">82dBm</a></literal></sticky>| <sticky><literal><a href="javascript:set_diag_high_txgain('-80');">80dBm</a></literal></sticky> | <sticky><literal><a href="javascript:set_diag_high_txgain('-78');">78dBm</a></literal></sticky>| <sticky><literal><a href="javascript:set_diag_high_txgain('-76');">76dBm</a></literal></sticky>| <sticky><literal><a href="javascript:set_diag_high_txgain('-74');">74dBm</a></literal></sticky>|


|<sticky><literal> <a  href="javascript:location.href=summary_vna_high_diag;"> <img height=480  src="" name="summary_vna_high_diag"></a> </literal></sticky>  |


<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_high_diagnostic_plots();
-->
</script>
</literal></sticky>
 
---++ Q vs T
Taking measurements of Q vs T below 1K. P = -100dBm. 

Logfile: ~/data/20180413_netanal/autonetanal_20180415_164043.txt

| <img src="%ATTACHURLPATH%/reso_86_to_562mK_QivsT.png" alt="reso_86_to_562mK_QivsT.png" width="500"  /> | <img src="%ATTACHURLPATH%/reso_86_to_562mK_dfvsT.png" alt="reso_86_to_562mK_dfvsT.png" width="500"  /> | 

---+++ VNA scans
<literal>
<script type='text/javascript'>// <![CDATA[
temp_frequency = '253.3';
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/tempsweep_plots/';

function update_temp_vna_plots(){
    summary_temp_vna_link = full_dir+ "rawfig/reso_"+temp_frequency+"MHz.png";

    document["summary_temp_vna"].src=summary_temp_vna_link;
}

function set_temp_frequency(x){
    temp_frequency = x;
    update_temp_vna_plots();
}
// ]]></script>
</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_temp_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_frequency('394.3');">394.3MHz</a></literal></sticky>|

| <sticky><literal> <a  href="javascript:location.href=summary_temp_vna;"> <img height=480  src="" name="summary_temp_vna"></a> </literal></sticky>  |

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_temp_vna_plots();
-->
</script>
</literal></sticky>

---+++ diagnostic fit plots
<literal>
<script type='text/javascript'>// <![CDATA[
temp_diag_frequency = '253.3';
temp_diag_temp = '86.1'
 
full_dir = '/twiki/pub/Resonators/20180413_TLS_Screener_NetAnal/tempsweep_plots/';

function update_temp_vna_diagnostic_plots(){
    summary_temp_vna_diag_link = full_dir+ "rawfig/diagnostic/reso_"+temp_diag_frequency+"MHz_"+temp_diag_temp+"mK.png";

    document["summary_temp_vna_diag"].src=summary_temp_vna_diag_link;
}

function set_temp_diag_temp_frequency(x){
    temp_diag_frequency = x;
    update_temp_vna_diagnostic_plots();
}

function set_temp_diag_temp(x){
    temp_diag_temp = x;
    update_temp_vna_diagnostic_plots();
}
// ]]></script>
</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('253.3');">253.3MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('254.2');">254.2MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('254.5');">254.5MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('283.7');">283.7MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('283.9');">283.9MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('286.0');">286.0MHz</a></literal></sticky>|<sticky><literal><a href="javascript:set_temp_diag_temp_frequency('394.3');">394.3MHz</a></literal></sticky>|

| <sticky><literal><a href="javascript:set_temp_diag_temp('86.1');">86.1mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('87.8');">87.8mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('96.7');">96.7mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('115.0');">115.0mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('138.4');">138.4mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('162.6');">162.6mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('188.6');">188.6mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('215.7');">215.7mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('242.1');">242.1mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('269.6');">269.6mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('297.8');">297.8mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('329.3');">329.3mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('357.7');">357.7mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('386.3');">386.3mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('414.9');">414.9mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('443.3');">443.3mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('472.6');">472.6mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('501.4');">501.4mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('531.0');">531.0mK</a></literal></sticky> | <sticky><literal><a href="javascript:set_temp_diag_temp('562.4');">562.4mK</a></literal></sticky> |

| <sticky><literal> <a  href="javascript:location.href=summary_temp_vna_diag;"> <img height=480  src="" name="summary_temp_vna_diag"></a> </literal></sticky>  |

<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_temp_vna_diagnostic_plots();
-->
</script>
</literal></sticky>

-- %USERSIG{AlbertWandui - 2018-05-18}%


