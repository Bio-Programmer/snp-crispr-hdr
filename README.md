# snp-crispr-hdr
<p>Make sure you include full file paths for all scripts below.</p>
snp-crispr-hdr.py
---------------
<p>Given a list of SNPs, determines which pass the CRISPR system constraints for HDR (homology-directed repair).</p>
<p>Command to run:</p>
<pre><code>python3 snp-crispr-hdr.py your_mpra_input.csv your_hdr_output.csv</pre></code>
crispr-motifbreakR.py 
---------------
<p>Given the output from the snp-crispr-hdr.py script and output from motifbreakR script (to be finalized), annotate your snps with motifbreakR data.</p>
<p>Command to run:</p>
<pre><code>python3 crispr-motifbreakR.py your_motifbreakR_input.csv your_snp-crispr-hdr_input your_crispr-motifbreakR_output.csv</pre></code>
TODO:
---------------
<p>Upload a clean motifbreakR script</p>
<p>Later on: enable this pipeline to generate the motifbreakR file without the user having to generate their own.</p>
