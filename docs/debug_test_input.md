### A to C ($\geq$ 1 fragment/complex)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type left --numfrag 1 --anchor 0 --output_file out_1
```

### A to C ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type left --numfrag 2 --anchor 0 --output_file out_2
```

### C to A ($\geq$ 1 fragment/complex)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type right --numfrag 1 --anchor 0 --output_file out_3
```

### C to A ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type right --numfrag 2 --anchor 0 --output_file out_4
```

### B to C ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type left --numfrag 2 --anchor 1 --output_file out_5
```

### B to A ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type right --numfrag 2 --anchor 2 --output_file out_6
```

### B towards A & C ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --path2 test_input.domains --type middle --numfrag 2 --anchor 3 --region chr3:150000-155000 --output_file out_7
```

### A and C and D ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --type multiple --numfrag 2 --region chr3:100000-108000\;chr3:300000-308000\;chr3:420000-428000 --operation yes\;yes\;yes --output_file out_8
```

### A, not B, and C and D ($\geq$ 2 fragments/complexes)
```
python main.py --path1 test_input.PEanno --type multiple --numfrag 2 --region chr3:100000-108000\;chr3:150000-155000\;chr3:300000-308000\;chr3:420000-428000 --operation yes\;no\;yes\;yes --output_file out_9
```

### Medium testing
```
python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type left --numfrag 2 --anchor 4 --output_file out_test_left_1

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type right --numfrag 2 --anchor 4 --output_file out_test_right_1

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --type multiple --numfrag 2 --region chr1:980316-988316\;chr1:1059964-1067964 --operation yes\;yes --output_file out_test_both_1

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type middle --numfrag 2 --anchor 4 --region chr1:1039000-1040000 --output_file out_test_middle_1
```

```
python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type middle --numfrag 2 --anchor 5 --region chr2:37654000-37657000 --output_file b_centered_test

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type left --numfrag 2 --anchor 5 --output_file left_test

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type right --numfrag 2 --anchor 5 --output_file right_test

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --region chr2:37392070-37400070\;chr2:37887111-37895111 --operation yes\;yes --type multiple --numfrag 2 --output_file both_test
```