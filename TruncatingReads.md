<!DOCTYPE html>
<html>
  <head>
      <meta charset="utf-8" />
      
  </head>
  <body class='markdown-preview' data-use-github-style><h1>Improving mergepairs with <a href="https://www.drive5.com/usearch/manual/cmd_fastx_truncate.html">fastx_truncate</a></h1>
<p>If you have poor merge success (&lt;60% merge), you may want to trim the reverse reads. Reverse reads are often of lower quality. Usearch merges by taking bp from the forward read, so removing bad basecalls in the reverse read can improve merging.</p>
<h4 id="you-can-use-the-mergepairs-alnout-file-to-determine-the-the-right-length-of-trimming-">you can use the mergepairs -alnout file to determine the the right length of trimming.</h4>
<p>For example, we trimmed the tail of the reverse read by 80bp because the forward read alignment to the reverse read tail was poor and the 505F to 815R has a higher overlap (i.e. ~210bp overlap)</p>
<h4 id="use-a-shell-for-loop-to-truncate">Use a shell for loop to truncate</h4>
<h4>You must be in the directory with you raw reverse reads</h4>
<pre class="editor-colors lang-"><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>mkdir&nbsp;trim_reverse_run1</span></span></div><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>nano&nbsp;rev_trunc_run1.sh</span></span></div></pre><pre class="editor-colors lang-"><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>for&nbsp;fq&nbsp;in&nbsp;*R2_001.fastq</span></span></div><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>do</span></span></div><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>&nbsp;&nbsp;/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64&nbsp;-fastx_truncate&nbsp;$fq&nbsp;-stripright&nbsp;80&nbsp;-fastqout&nbsp;trim_reverse_run1/$fq</span></span></div><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>&nbsp;</span></span></div><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>done</span></span></div></pre><p>Make the file executable and run the code.</p>
<pre class="editor-colors lang-"><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>chmod&nbsp;+x&nbsp;nano&nbsp;rev_trunc_run1.sh</span></span></div><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>./nano&nbsp;rev_trunc_run1.sh</span></span></div></pre><p>Now try to merge the sequences again with the truncated reverse reads.</p>
<h4>note that you need to specify the location for your forward (R1) and reverse reads (R2)</h4>
<pre class="editor-colors lang-"><div class="line"><span class="syntax--text syntax--plain syntax--null-grammar"><span>/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64&nbsp;-fastq_mergepairs&nbsp;*R1*.fastq&nbsp;-reverse&nbsp;rim_reverse_run1/*R2*&nbsp;-relabel&nbsp;@&nbsp;-fastqout&nbsp;mergedfastq_run1/combined_merged_run1_rev_trun.fastq&nbsp;&nbsp;-tabbedout&nbsp;mergedfastq_run1/combined_merged_run1_rev_trun_pair_report.txt&nbsp;-alnout&nbsp;mergedfastq_run1/combined_merged_run1_rev_trun_pair_aln.txt</span></span></div></pre><p>The truncating of the reverse reads by <strong>80bp</strong> increased our merging from <strong>33.29%</strong> to <strong>67.2%</strong></p>
<p>If this does not improve the percentage merged, consider just using the forward reads as recommended by <a href="https://drive5.com/usearch/manual/merge_badrev.html">Edgar</a></p></body>
</html>
