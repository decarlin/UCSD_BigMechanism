targets: diffused.txt

stem_cell_2: MEP GMP
	rm *.tmp

MEP:
	../automateHeatKernel.py -u decarlin -n aa10eb78-8278-11e5-b435-06603eb7f303 -p perfect6 -q changed_MEP.txt -x -s MEP_changed.sif -o MEP_filtered.sif -d diffused_MEP.tmp -N 100 -k kernel.tmp

GMP:
	../automateHeatKernel.py -u decarlin -n aa10eb78-8278-11e5-b435-06603eb7f303 -p perfect6 -q changed_GMP.txt -x -s GMP_changed.sif -o GMP_filtered.sif -d diffused_GMP.tmp -N 100 -k kernel.tmp

erna_phospho:
	../automateHeatKernel.py -q Erna_phospho_perturbed.txt -H ../mappings/bel_HGNCtoRGD.tab -m ../mappings/bel_MGItoRGD.tab -s erna_phospho.sif -d diffused_neural_repair_phospho.txt -N 200 -S 'rat'

established_kernel: some_genes.txt
	../multiplyEstablishedKernel.py -k kernel.txt -q some_genes.txt -d diffused_established.txt

stem_example_human.txt: human_GO_regulation_of_stem_cell_maintenence.txt
	../automateHeatKernel.py -q human_GO_stem_cell_maintenence.txt -m ../mappings/bel_MGItoHGNC.tab -r ../mappings/bel_RGDtoHGNC.tab -s stem_example_human.sif -d diffused_stem_human.txt -S 'human'

mouse_neural_repair_example.txt: mouse_neural_repair.txt
	../automateHeatKernel.py -q mouse_neural_repair.txt -H ../mappings/bel_HGNCtoMGI.tab -r ../mappings/bel_RGDtoMGI.tab -s stem_example.sif -d diffused_stem.txt -S 'mouse'

pancrease_wighted_example:
	cut.pl -f 1,3 ~/projects/cell_context_mask/GTEx_pancreas_rpkm_stats.tab \
	| join.pl -o 'NA' diffused_pancreas_weighted.txt - \
	> pancrease_wighted_example.tab

diffused.txt: some_genes2.txt
	../automateHeatKernel.py -q some_genes2.txt -m ../mappings/bel_MGItoHGNC.tab -r ../mappings/bel_RGDtoHGNC.tab

nonbel_example: some_genes.txt
	../automateHeatKernel.py -q some_genes.txt -m ../mappings/bel_MGItoHGNC.tab -r ../mappings/bel_RGDtoHGNC.tab -n d106d033-baf0-11e5-abf0-06603eb7f303 -s nonbel.sif -k nonbel_kernel.txt

weighted_query: some_weighted_genes.txt
	../automateHeatKernel.py -w some_weighted_genes.txt -m ../mappings/bel_MGItoHGNC.tab -r ../mappings/bel_RGDtoHGNC.tab -b weighted_out.bel

several_network_query: some_big_networks.txt
	../multipleHeatKernel.py -q some_genes.txt -m ../mappings/bel_MGItoHGNC.tab -r ../mappings/bel_RGDtoHGNC.tab -s multiple_nets.sif -n some_big_networks.txt -c kernel

full_kernel:
	../buildFullKernel.py -k ../kernels/BEL_Large_Corpus_HGNC.txt -m ../mappings/bel_MGItoHGNC.tab -r ../mappings/bel_RGDtoHGNC.tab