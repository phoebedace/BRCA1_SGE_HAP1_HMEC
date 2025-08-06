#run_move_NS_from_asf.py
import os
import sys
import subprocess

def xlsx(fname):
    import zipfile
    from xml.etree.ElementTree import iterparse
    z = zipfile.ZipFile(fname)
    strings = [el.text for e, el in iterparse(z.open('xl/sharedStrings.xml')) if el.tag.endswith('}t')]
    rows = []
    row = {}
    value = ''
    for e, el in iterparse(z.open('xl/worksheets/sheet1.xml')):
        if el.tag.endswith('}v'):  # Example: <v>84</v>                            
            value = el.text
        if el.tag.endswith('}c'):  # Example: <c r="A3" t="s"><v>84</v></c>                                 
            if el.attrib.get('t') == 's':
                value = strings[int(value)]
            letter = el.attrib['r']  # Example: AZ22                         
            while letter[-1].isdigit():
                letter = letter[:-1]
            row[letter] = value
            value = ''
        if el.tag.endswith('}row'):
            rows.append(row)
            row = {}
    return rows

#usage python ~/home/users/findlag/bin/run_move_from_asf.py $seq_dir
#out_dir must contain the .xlsx file with sample names provided by asf
#will output a script to be run in shell

working_dir = os.getcwd()
shell_file = open('run_move_from_asf.sh', 'w')
shell_file.write('cd '+working_dir+'\n') #makes sure shell script runs on correct folder only
seq_dir = sys.argv[1]

sample_sheets_found = 0
for i in os.listdir(seq_dir):
    if i.endswith('.xlsx'):
        sample_sheets_found+=1
        sample_sheet_data = xlsx(seq_dir+'/'+i)
    else:
        sample_sheet_data = xlsx(sys.argv[2])

if sample_sheets_found == 1:
    print("Sample sheeet found in run folder")
else:
    print("Requiring sample sheet to be added to output directory...")
    print("Sample sheets found in run directory: ")
    print(sample_sheets_found)

#needs to be written to open the .xlsx file in the run directory and retrieve sample pairings from it:

#print sample_sheet_data; sample_sheet imported in format list, with each row a dictionary pairing column (A,B,C, etc) to entry
#Important formatting: will skip two lines, and then make a pairing for all entries based on A and B (first two columns)

sample_dict = {}
lines_counted = 0
for xlsx_line in sample_sheet_data:
    lines_counted +=1
    if lines_counted < 2:
        pass
    else:
        sample_dict[xlsx_line['A']] = xlsx_line['B']

#now have list of pairings to rename

#reference comes first, and then the cigar reflects changes from reference to sample
for i in os.listdir(seq_dir+'/fastq'):
    if i.endswith("fastq.gz"): #will only work on fastq.gz files
        index_of_first_S = i.find('_S')
        index_of_first_under_post_S = i.find('_', index_of_first_S+1)
        old_name = i[:index_of_first_S]
        if old_name in sample_dict:
            sample_name = sample_dict[old_name]+i[index_of_first_under_post_S:]
            shell_file.write("cp "+seq_dir+'/fastq/'+i+' fastq/'+sample_name+'\n')
        else:
            pass
    else:
        pass

for i in os.listdir(seq_dir+'/fastqc'):
    if i.endswith("fastqc.html"): #will only work on fastqc.htmlfiles
        index_of_first_S = i.find('_S')
        index_of_first_under_post_S = i.find('_', index_of_first_S+1)
        old_name = i[:index_of_first_S]
        if old_name in sample_dict:
            sample_name = sample_dict[old_name]+i[index_of_first_under_post_S:]
            shell_file.write('cp '+seq_dir+'/fastqc/'+i+' fastqc/'+sample_name+'\n')
            #example ending FIN2752A25_S18_L004_R2_001_fastqc.html
        else:
            pass
    else:
        pass

shell_file.close()

#script format:  python