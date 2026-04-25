input_file = "/home/charlie/Downloads/mm39_NCBI_RefSeq.fa"
output_file = open("/home/charlie/Downloads/mm39_NCBI_RefSeq_no_dups.fa", 'w')
removed_lines = open('removed_lines.fa', 'w')
names = set({})
printing = True
for line in open(input_file, 'r'):
	if ">" in line:
		line = line.split(" ")[0].replace("mm39_ncbiRefSeqCurated_","")+"\n"
		if line not in names:
			printing = True
		if line in names:
			printing = False
		names.add(line)
	if printing == True:
		output_file.write(line)
	if printing == False:
		removed_lines.write(line)

