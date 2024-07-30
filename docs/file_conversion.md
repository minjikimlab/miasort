## pairs2regions

**7m46.145s**

58,649,978 lines

```Shell
time ./pairs2regions . ./LHG0035N_0035V_0045V.bsorted.pairs.gz ./hg38.chrom.sizes
```

with no batch processing nor threadpool, 7 mins 46s
with no batch processing, 6 mins

8-Core Intel I9

batch_size = 100
96s

batch_size = 1,000
79s

batch_size = 10,000
77s

batch_size = 100,000
75s

batch_size = 1,000,000
81s

````
(env) zhangzzc@mac-c02yk0jqjv3y complex-sorter % ./bin/time-testing-hic.sh
+ g++ -std=c++11 -o pairs2regions pairs2regions.cpp -lz
++ date +%s
+ start_time=1722358678
+ ./pairs2regions . ./4DNFIACOTIGL.pairs.gz ./hg38.chrom.sizes

real    24m6.579s
user    161m12.569s
sys     4m21.449s
++ date +%s
+ end_time=1722360125
+ duration=1447
+ echo 'Total time for running pairs2regions: 1447 seconds'
Total time for running pairs2regions: 1447 seconds
```