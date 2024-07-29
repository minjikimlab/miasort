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