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