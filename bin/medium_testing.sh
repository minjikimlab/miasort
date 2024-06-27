python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type left --numfrag 2 --anchor 4 --output_file out_test_left_1

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type right --numfrag 2 --anchor 4 --output_file out_test_right_1

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --type multiple --numfrag 2 --region chr1:980316-988316\;chr1:1059964-1067964 --operation yes\;yes --output_file out_test_both_1

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test_input.domains --type middle --numfrag 2 --anchor 4 --region chr1:1039000-1040000 --output_file out_test_middle_1