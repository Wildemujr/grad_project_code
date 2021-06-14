#!/usr/bin/env python3
 
from pprint import pprint as pp
import argparse
import numpy as np
import pandas as pd
import requests



def get_kegg_ko_id(kegg_gene_ids):
	
	orthology_ids = [] 
	for gene_id in kegg_gene_ids:
		r = requests.get(f"http://rest.kegg.jp/get/{gene_id}")
		content = r.text
		content_list = content.split("\n")
		for line in content_list:
			if "ORTHOLOGY" in line:
				print(line)
				temp_l = line.split(" ")
				temp_l = list(filter(None, temp_l))
				orthology_ids.append((gene_id, temp_l[1]))
				break

	return orthology_ids




if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-i",
		"--input",
		dest="inputf",
		help="Choosing the Input File to be Processed and Reformatted."
	)
	parser.add_argument(
		"-o",
		"--output",
		dest="outputf",
		help="Choosing the Name of the Ouput File to be Output After Processing is Complete."
	)
	option = parser.parse_args()
	
	# Reading in Excel sheet as a dataframe and promptly converting to list.
	df = pd.read_excel(f"{option.inputf}").iloc[:,6].values.tolist()
	
	# Converting all list elements to human kegg gene id's
	kegg_ids = list(map(lambda x:f"hsa:{x}",df))
	ids_comb = list(zip(df, kegg_ids))
	kegg_ko_ids = get_kegg_ko_id(kegg_ids)
	# print(kegg_ids)	

	f = open("kegg_overlap_input.txt","w+")
	for i in kegg_ko_ids:
		f.write(f"{i[0]}  {i[1]}\n")
	f.close()

	entrez_kegg_df = pd.DataFrame(
		data=ids_comb,
		columns=["Entrez ID","KEGG ID"]
	)
	entrez_kegg_df.to_excel(f"{option.outputf}.xlsx")