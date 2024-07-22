<div align="center">
  <h1 align="center">ComplexSorter</h1>
</div>

## Installation
```Shell
chmod +x ./bin/install.sh
./bin/install.sh
source env/bin/activate
```

Bedtools should also be installed. To do so, follow instructions on https://bedtools.readthedocs.io/en/latest/content/installation.html#.
One of the options is:
```Shell
$ mamba install -c conda-forge bedtools
```

## Example Command for Running the Tool

```Shell
python main.py \
--path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno \
--path2 test-july-2.bedte --type abc --graphs BtoA --numfrag_min 2 \
--anchor_options no --out_dir test_folder
```

## Input Files
### Anchors and Regions (`.domains`)

```
chr1	982867	987410	chr1	1371055	1372952	1	313	363	984501	1372412	dm0	chr1	984316	984335	+	185	chr1	1372224	1372243	-	188
chr1	2194592	2197930	chr1	2378563	2382994	1	602	316	2195215	2381889	dm1	chr1	2195348	2195367	+	133	chr1	2381915	2381934	-	26
chr1	3452263	3454407	chr1	3618470	3620776	1	529	715	3453278	3619242	dm2	chr1	3453290	3453309	+	12	chr1	3619183	3619202	-	59
```

Every line is independent from each other. We will only do one line at a time. Each line will be conducted in only 1 second.

Column 1 and 2 (0-indexed) indicate the range of **left anchor** (left reference) and Column 4 and 5 indicate the **right anchor**.

We will need right anchor even when we do left anchor because we only want to rank GEMs that meet the left anchor and also with the rightmost fragment that is less than the right anchor region.

### GEMs (`.PEanno`)

```
chr10	101762606	101763188	2	SHG8033-100000080-AAACACCAGAGGGTAABX8033-HEA-5-2-sub-1-1	E
chr10	104102286	104102852	2	SHG8033-100000080-AAACACCAGAGGGTAABX8033-HEA-5-2-sub-1-1	E
chr10	88838102	88838730	2	SHG8033-100000171-AAACACCAGCGCCCTABX8033-HEA-4-0-sub-1-1	P
```

**Every line represents a fragment**. Column 1 (0-indexed) represents the left point and Column 2 represents the right point of the fragment.

Column 3 (0-indexed), in this case `2`, indicates **the number of fragments in one GEM/chromatin complex**.

**Each GEM has a unique id** (Column 4), so you can see that line 1 and 2 have the same Column 4 since they are the two fragments in the same GEM.

Ignore the 6th column.

## License
Shield: [![CC BY-NC-ND 4.0][cc-by-nc-nd-shield]][cc-by-nc-nd]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International License][cc-by-nc-nd].

[![CC BY-NC-ND 4.0][cc-by-nc-nd-image]][cc-by-nc-nd]

[cc-by-nc-nd]: http://creativecommons.org/licenses/by-nc-nd/4.0/
[cc-by-nc-nd-image]: https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png
[cc-by-nc-nd-shield]: https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg
