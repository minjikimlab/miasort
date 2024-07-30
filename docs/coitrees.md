```
(env) zhangzzc@mac-c02yk0jqjv3y complex-sorter % ./bin/coitrees-testing.sh
++ date +%s
+ start_time=1722270867
+ cargo run --release --example bed-intersect -- GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno test-july-2.bedte
    Finished `release` profile [optimized] target(s) in 0.10s
     Running `target/release/examples/bed-intersect GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno test-july-2.bedte`
reading bed: 4.436s
lines: 42843540
sequences: 24
veb_order: 4.227s
overlap: 0.004s
total overlaps: 569

real    0m9.330s
user    0m6.166s
sys     0m1.252s
++ date +%s
+ end_time=1722270876
+ duration=9
+ echo 'Total time for running coitrees intersect: 9 seconds'
Total time for running coitrees intersect: 9 seconds
```

```
(env) zhangzzc@mac-c02yk0jqjv3y complex-sorter % ./bin/coitrees-testing.sh
++ date +%s
+ start_time=1722360225
+ cargo run --release --example bed-intersect -- 4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region test-july-2.bedte
   Compiling coitrees v0.4.0 (/Users/zhangzzc/Desktop/complex-sorter)
    Finished `release` profile [optimized] target(s) in 21.10s
     Running `target/release/examples/bed-intersect 4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region test-july-2.bedte`
reading bed: 166.754s
lines: 317454742
sequences: 24
veb_order: 35.903s
overlap: 0.024s
total overlaps: 17

real    3m45.496s
user    0m47.274s
sys     0m17.035s
++ date +%s
+ end_time=1722360451
+ duration=226
+ echo 'Total time for running coitrees intersect: 226 seconds'
Total time for running coitrees intersect: 226 seconds
```