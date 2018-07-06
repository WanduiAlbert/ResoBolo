<literal>
<script type='text/javascript'>// <![CDATA[
IQ_frequency = '320';
IQ_txgain = '0';

full_dir = '/twiki/pub/Resonators/20180703_TLS_Screener_NEF/spectra/';

function update_vna_IQ_noise_comp_plots(){
    summary_IQ_vna_link = full_dir+ "centered_IQ_"+IQ_frequency+"MHz_txgain_"+ IQ_txgain+ "dB.png";

    document["summary_IQ_vna_comp"].src=summary_IQ_vna_link;
}

function set_IQ_frequency(x){
    IQ_frequency = x;
    update_vna_IQ_noise_comp_plots();
}

function set_IQ_txgain(x){
    IQ_txgain = x;
    update_vna_IQ_noise_comp_plots();
}
// ]]></script>
</literal> <sticky><table width="100%"><tr><td valign="top" width="25%"></sticky>

| <sticky><literal><a href="javascript:set_IQ_frequency('320');">320MHz</a></literal></sticky> | <sticky><literal><a href="javascript:set_IQ_frequency('360');">360MHz</a></literal></sticky> |

| <sticky><literal><a href="javascript:set_IQ_txgain('0');">0 dB</a></literal></sticky> | <sticky><literal><a href="javascript:set_IQ_txgain('5');">5 dB</a></literal></sticky> | <sticky><literal><a href="javascript:set_IQ_txgain('10');">10 dB</a></literal></sticky> | <sticky><literal><a href="javascript:set_IQ_txgain('15');">15 dB</a></literal></sticky> | <sticky><literal><a href="javascript:set_IQ_txgain('20');">20 dB</a></literal></sticky> | <sticky><literal><a href="javascript:set_IQ_txgain('25');">25 dB</a></literal></sticky> |

| <sticky><literal> <a  href="javascript:location.href=summary_IQ_vna_comp;"> <img height=600  src="" name="summary_IQ_vna_comp"></a> </literal></sticky>  |
<sticky></td></tr></table>

<literal>
<script language="JavaScript">
<!--
update_vna_IQ_noise_comp_plots();
-->
</script>
</literal></sticky>  