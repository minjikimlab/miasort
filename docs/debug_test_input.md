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