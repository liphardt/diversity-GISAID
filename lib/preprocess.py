import pandas as pd
from datetime import date
from datetime import timedelta
import argparse
import warnings
import os
from Bio import SeqIO
warnings.simplefilter(action='ignore', category=FutureWarning)


#function to fix dates in metadata to proper iso format
def fix_time(input_frame):
	input = input_frame.to_dict('records')
	output_frame = pd.DataFrame(columns = input_frame.columns)
	for row in input:
		if len(row['date'].split("-")) < 3:
			continue
		if len(row['date'].split("-")[1]) < 2 and len(row['date'].split("-")[2]) < 2:
			day = f"0{row['date'].split('-')[2]}"
			month = f"0{row['date'].split('-')[1]}"
			year = row['date'].split("-")[0]
			date_list = [year, month, day]
			new_date = '-'.join(date_list)
			row['date'] = new_date
			output_frame = output_frame.append(pd.DataFrame(row, index = [1]),\
 			ignore_index = True)
		if len(row['date'].split("-")[2]) < 2:
			day = f"0{row['date'].split('-')[2]}"
			year_month = row['date'].split("-")[:2]
			year_month.append(day)
			new_date = '-'.join(year_month)
			row['date'] = new_date
			output_frame = output_frame.append(pd.DataFrame(row, index = [1]),\
 			ignore_index = True)
		if len(row['date'].split("-")[1].strip()) < 2:
			month = f"0{row['date'].split('-')[1]}"
			year = row['date'].split("-")[0]
			day = row['date'].split("-")[2]
			date_list = [year, month, day]
			new_date = '-'.join(date_list)
			row['date'] = new_date
			output_frame = output_frame.append(pd.DataFrame(row, index = [1]),\
 			ignore_index = True)
		else:
			output_frame = output_frame.append(pd.DataFrame(row, index = [1]),\
 			ignore_index = True)
	return output_frame

#Add relevant data fields we want to populate. 
def add_data_columns(input_frame):
	input_frame['time_period'] = ''
	input_frame['with_locality'] = ''
	return input_frame

#Count total sequences for weekly time periods across all labs, add time period for each seq
def count_total(input_frame):
	week = timedelta(days=7)
	meta_dates = sorted(input_frame.date)
	temp_frame = pd.DataFrame(columns = ["period", "count_total", "count_locality"])
	begin = date.fromisoformat("2020-01-01")
	ending = date.fromisoformat(meta_dates[-1])
	while begin < ending:
		period = begin + week
		count_total = 0
		count_locality = 0
		for row in input_frame.itertuples():
			if begin <= date.fromisoformat(row.date) <= period:
				index = input_frame.loc[input_frame.strain == row.strain].index
				input_frame.at[index[0], 'time_period'] = f"{begin}_{period}"
				count_total += 1
				if pd.notnull(row.location):
					count_locality += 1
		temp_frame = temp_frame.append({'period': f"{begin}_{period}",
			"count_total": count_total, "count_locality": count_locality}, ignore_index=True)
		begin += week
	file_name = "results/all_sample_period_counts.csv"
	temp_frame.to_csv(file_name, index = False)

#Add whether seq has locality data or not
def check_loc(input_frame):
	no_count = 0	
	for row in input_frame.itertuples():
		index = input_frame.loc[input_frame.strain == row.strain].index
		if pd.isnull(row.location):
			input_frame.at[index[0], 'with_locality'] = 'no'
			no_count += 1
		else:
			input_frame.at[index[0], 'with_locality'] = 'yes'
	return input_frame

#Count seq for each period for each lab of and write new csv
def count_period_cases(lab, case_data):
	meta_dates = sorted(metadata.date)
	temp_frame = pd.DataFrame(columns = ["period", "count"])
	temp_dict = case_data.to_dict('records')
	begin = date.fromisoformat(meta_dates[1])
	ending = date.fromisoformat(meta_dates[-1])
	while begin < ending:
		period = begin + week
		count = 0
		for row in temp_dict:
			if row['submitting_lab'] == lab:
				if begin <= date.fromisoformat(row['date']) <= period:
					count += 1
		temp_frame = temp_frame.append({'period': f"{begin}_{period}", "count": count}, ignore_index=True)
		begin += week
	file_name = lab.split(" ")[0] + ".csv"
	temp_frame.to_csv(file_name, index = False)

#Add WHO clades to spreadsheet
def add_who(input_frame):
	who_clades = {"GRY":"Alpha", "GH": "Beta", "GR":"Gamma", "GK": "Delta",\
				"GRA": "Omicron"}
	clade_list = []
	for row in input_frame.itertuples():
		if row.GISAID_clade in who_clades:
			clade_list.append(who_clades[row.GISAID_clade])
		else:
			clade_list.append("Other")
	position = len(input_frame.columns) - 1
	input_frame.insert(position, "who_clades", clade_list)
	return input_frame


#Calculate for each lab	the number of sequences with and without location data. 
def count_labs(input_frame):
	temp_dict = input_frame.to_dict('records')
	sub_lab = set(input_frame.submitting_lab)
	temp_frame = pd.DataFrame(columns = ['lab', 'count_total', 'count_locality'])
	for lab in sub_lab:
		count_total = 0
		count_locality = 0
		for row in temp_dict:
			if row['submitting_lab'] == lab:
				count_total += 1
				if row['with_locality'] == 'yes':
					count_locality += 1
		temp_frame = temp_frame.append(pd.DataFrame({'lab': lab, 'count_total': count_total, 'count_locality': count_locality}, index = [1]), ignore_index = True)
	file_name = 'submitting_lab_counts.csv'
	temp_frame.to_csv(file_name, index = False)	

def pull_meta_seqs(meta_frame, sequences):
        fastas = list(SeqIO.parse(sequences, "fasta"))
        fasta_dict = {}
        for i in range(0,len(fastas)):
            fasta_dict[fastas[i].id] = fastas[i]
        metadict = meta_frame.to_dict('records')
        unique_dates = []
        parent_dir = "results/week_subset"
        for row in meta_frame.itertuples():
            unique_dates.append(str(row.time_period))
        unique_dates = set(unique_dates)
        for period in unique_dates:
            if period != "nan" and period != "":
                path = os.path.join(parent_dir, period)
                os.makedirs(path)
        for period in unique_dates:
            if period != "nan" and period != "":
                begin = date.fromisoformat(period.split("_")[0])
                end = date.fromisoformat(period.split("_")[1])
                filename = f"{parent_dir}/{period}/{period}_meta.csv"
                temp_meta = pd.DataFrame(columns = meta_frame.columns)
                for record in metadict:
                    if record['date'] != "nan":
                        if begin <= date.fromisoformat(record['date']) <= end:
                            temp_frame = pd.DataFrame(record, index = [0])
                            temp_meta = pd.concat([temp_meta, temp_frame], ignore_index = True)
                temp_meta.to_csv(filename, index = False)
        for root,dirs,files in os.walk(parent_dir):
            for file in files:
                full_file = os.path.join(root, file)
                month_data = pd.read_csv(full_file)
                temp_fasta = []
                period = file.split(".")[0]
                file_name = f"{period}.sequences.fasta"
                pathway = os.path.join(root,file_name)
                for row in month_data.itertuples():
                    seq = fasta_dict[row.strain]
                    if len(seq.seq) > 28000:
                        temp_fasta.append(seq)
                if len(temp_fasta) > 2:
                    with open(pathway, "w") as handle:
                        SeqIO.write(temp_fasta, handle, "fasta")

if __name__ == "__main__":
	#Set globals and read input data
	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_meta', type = str, help = 'Input file name for meta data from gisaid')
	parser.add_argument('--output_meta', type = str, help = 'output file name for processed meta data')
	parser.add_argument('--loc_data', type = str, help = 'Locality data for UMGC samples')
	parser.add_argument('--gen_labs', type = str, choices = ['yes'], help = 'Generate per lab CSVs of sequences per weekly period')
	parser.add_argument('--add_who_clades', type = str, choices = ['yes'], help = 'Add WHO clades of interest (like omicron, delta, etc.) from GISAID clade')
	args = parser.parse_args()
	metadata = pd.read_csv(args.input_meta)#, sep = "\t")
	metadata.sort_values(by = ['date'], inplace = True)
	outmeta = args.output_meta
	if args.fix_locality == 'yes':
		umgc_locality_data = pd.read_csv(args.loc_data)
	week = timedelta(days=7)
	
	
	print("Mmk, let's do this")
	print(f"Input is {args.input_meta}")
	print("Fixing time so it is ISO formatted")
	fixed_metadata = fix_time(metadata)
	print(f"Adding data fields")
	add_data_columns(fixed_metadata)
	print("Setting weekly time periods and counting cases")
	count_total(fixed_metadata)
	print("Checking whether seqs have locality data")
	check_loc(fixed_metadata)
	print("Counting seqs/period/lab and writing to csv")
	if args.gen_labs == 'yes':
		for lab in set(fixed_metadata.submitting_lab):
			print(f"Working on {lab}")
			count_period_cases(lab, fixed_metadata)
	if args.add_who_clades == 'yes':
		add_who(fixed_metadata)
	print("Writing updated metadata to file")
	fixed_metadata.to_csv(outmeta, index = False)
	print("Writing lab count csv")
	count_labs(fixed_metadata)
