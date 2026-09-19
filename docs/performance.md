# Performance

The September 19 benchmark used ViroSync 1.0.1 with GenomeTools seed indexing,
shared repeat comparisons and development resources v1.1.0. Its source SHA-256 was
`a19fe5516183507df1c5d16cf794670baef0dfa769fd78bf7f135ee53a2ec3e3`.
ViroSync ran with one worker and 16 threads per input, at most two inputs at
once, and no optional analyses. Select a figure to open the full-size image.

## Synthetic boundary recovery

The synthetic set contains 60 loci and 3,365.2 kb of inserted sequence. Mean
best-call boundary recall was 0.777 for ViroSync, 0.596 for ViralRecall v3.1.0,
0.307 for ViralRecall v2, 0.056 for DetectEVE v1.4.0, and 0.066 for EEfinder
v1.1.1.

[![Synthetic boundary recovery and call burden](assets/performance/benchmark_fig2_syn2_detection_performance.png)](assets/performance/benchmark_fig2_syn2_detection_performance.png)

*Boundary recovery and call burden across 60 synthetic loci. Panel a credits
the call with the highest Jaccard index at each locus. Panel b compares missed
and outside-truth sequence. Panel c partitions calls at a 500 bp one-to-one
threshold. Missed-locus labels use recovery after the union of call fragments,
which can recover a locus without a credited one-to-one call. The test excludes
30 survey-derived loci whose source elements ViroSync found before this
comparison. Region detection has no true negative,
so the plot does not report specificity or accuracy.*

## Runtime and peak memory

Mean wall time across 30 SynEVEs-2 inputs was 100.5 seconds for ViroSync and
12.3 seconds for ViralRecall v3.1.0. Across five real-genome inputs, the means
were 600.0 and 1,705.5 seconds. Campaign concurrency differed, so these values
do not give a general speed ranking.

[![Runtime and peak resident memory](assets/performance/benchmark_figS1_runtime_memory_with_vr30.png)](assets/performance/benchmark_figS1_runtime_memory_with_vr30.png)

*Per-input wall time and peak resident memory for 30 SynEVEs-2 inputs and five
real-genome inputs. All tools used 16 threads or cores on the same 64-core
host, but campaign concurrency differed. All ViroSync timings come from the
September 19 campaign. Wall time includes contention and database caching.
Peak RSS is shown in decimal GB for the largest single
process, not the process tree. EEfinder is absent because its campaign used
different load and concurrency.*

## Real-genome output

ViroSync accepted 524 regions across the five real genomes. These inputs have
no complete region-level truth set. The figure describes output burden and
gene content, not detection accuracy.

[![Real-genome candidate burden and gene composition](assets/performance/benchmark_fig3_real_burden_composition.png)](assets/performance/benchmark_fig3_real_burden_composition.png)

*Candidate count, summed region length, and gene composition across five real
genomes. DetectEVE and EEfinder emit short homology hits. Their counts measure
output granularity and are not direct locus-count comparisons.*

See [Outputs](METHODS.md) for the result files and their interpretation.
