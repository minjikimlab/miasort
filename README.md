<div align="center">
  <h1 align="center">ComplexSort</h1>
</div>

## Run the Program
An example terminal command is:
```
python main.py
    --path1 CDH0002NR9L_hg38_CTCF_filt_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno
    --path2 first10.domains
    --type both
    --anchor 1
    --output_file out
```

## Input Files
### Anchors and Regions (`.domains`)

```
chr1	982867	987410	chr1	1371055	1372952	1	313	363	984501	1372412	dm0	chr1	984316	984335	+	185	chr1	1372224	1372243	-	188
chr1	2194592	2197930	chr1	2378563	2382994	1	602	316	2195215	2381889	dm1	chr1	2195348	2195367	+	133	chr1	2381915	2381934	-	26
chr1	3452263	3454407	chr1	3618470	3620776	1	529	715	3453278	3619242	dm2	chr1	3453290	3453309	+	12	chr1	3619183	3619202	-	59
chr1	5508750	5512637	chr1	5744829	5747502	1	272	263	5510287	5745950	dm3	chr1	5510051	5510070	+	236	chr1	5745972	5745991	-	22
chr1	6245977	6249005	chr1	6404310	6406063	1	418	328	6246688	6405008	dm4	chr1	6246886	6246905	+	198	chr1	6404836	6404855	-	172
chr1	6717812	6720523	chr1	7069291	7072256	1	455	459	6719940	7070555	dm5	chr1	6720164	6720183	+	224	chr1	7070576	7070595	-	21
chr1	7663117	7670765	chr1	7748190	7756109	1	380	733	7667816	7752892	dm6	chr1	7667705	7667724	+	111	chr1	7752995	7753014	-	103
chr1	8013338	8020511	chr1	8310262	8317605	1	890	334	8015209	8314479	dm7	chr1	8015170	8015189	+	39	chr1	8314472	8314491	-	7
chr1	9106882	9113411	chr1	9193303	9199911	1	642	284	9110922	9197002	dm8	chr1	9110927	9110946	+	5	chr1	9196621	9196640	-	381
chr1	9266080	9270035	chr1	9489154	9495246	1	484	530	9267276	9491249	dm9	chr1	9267277	9267296	+	1	chr1	9491238	9491257	-	11
```

Every line is independent from each other. We will only do one line at a time. Each line will be conducted in only 1 second.

Column 1 and 2 (0-indexed) indicate the range of **left anchor** (left reference) and Column 4 and 5 indicate the **right anchor**.

We will need right anchor even when we do left anchor because we only want to rank GEMs that meet the left anchor and also with the rightmost fragment that is less than the right anchor region.

### GEMs (`.PEanno`)

```
chr10	101762606	101763188	2	SHG8033-100000080-AAACACCAGAGGGTAABX8033-HEA-5-2-sub-1-1	E
chr10	104102286	104102852	2	SHG8033-100000080-AAACACCAGAGGGTAABX8033-HEA-5-2-sub-1-1	E
chr10	88838102	88838730	2	SHG8033-100000171-AAACACCAGCGCCCTABX8033-HEA-4-0-sub-1-1	P
chr10	88862713	88863341	2	SHG8033-100000171-AAACACCAGCGCCCTABX8033-HEA-4-0-sub-1-1	P
chr10	10867803	10868431	2	SHG8033-100000309-AAACACCAGGTCTCGCBX8033-HEA-6-3-sub-1-1	E
chr10	11012470	11013098	2	SHG8033-100000309-AAACACCAGGTCTCGCBX8033-HEA-6-3-sub-1-1	E
chr10	47464972	47465589	2	SHG8033-100000518-AAACACCCAATAGCGGBX8033-HEA-6-4-sub-1-1	E
chr10	47648402	47649030	2	SHG8033-100000518-AAACACCCAATAGCGGBX8033-HEA-6-4-sub-1-1	E
chr10	114112356	114112984	2	SHG8033-100000650-AAACACCCACTCACCTBX8033-HEA-4-0-sub-1-1	E
chr10	114217837	114218465	2	SHG8033-100000650-AAACACCCACTCACCTBX8033-HEA-4-0-sub-1-1	E
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