import sys
import csv

def read_lines_and_create_lists(filename):
    data_list = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        
        if len(lines) >= 11314:
            for i in range(76,11314):
                temp_line = lines[i].strip().split('\t')
                
                filename_without_extension = filename.split('.')[0]
            
                temp_line.insert(0, filename_without_extension)
            
                data_list.append(temp_line)
    
    return data_list

def process_list(data_list):
    new_list = [["SampleID", "position", "ref_base", "depth", "count_A", "count_C", "count_G", "count_T"]]
    
    for line in data_list:
        new_line = ["None"] * 8
        new_line[0] = line[0]
        new_line[1] = line[2]
        new_line[2] = line[3]
        new_line[3] = line[4]
        for element in line[6:]:
            element_parts = element.split(':')
            strand_bias = int(element_parts[5])/int(element_parts[1]) * 100 
            if element_parts[0] == 'A' and strand_bias >= 10 and strand_bias <= 90:
                #new_line[4] = str(round(int(element_parts[1]) / float(line[4]) * 100, 2))
                new_line[4] = element_parts[1]
            elif element_parts[0] == 'C'and strand_bias >= 10 and strand_bias <= 90:
                #new_line[5] = str(round(int(element_parts[1]) / float(line[4]) * 100, 2))
                new_line[5] = element_parts[1]
            elif element_parts[0] == 'G'and strand_bias >= 10 and strand_bias <= 90:
                #new_line[6] = str(round(int(element_parts[1]) / float(line[4]) * 100, 2))
                new_line[6] = element_parts[1]
            elif element_parts[0] == 'T'and strand_bias >= 10 and strand_bias <= 90:
                #new_line[7] = str(round(int(element_parts[1]) / float(line[4]) * 100, 2))
                new_line[7] = element_parts[1]
            else:
                this = "means nothing"

        new_list.append(new_line)
    
    return new_list

# Get input tsv as the first command line arg
filename = sys.argv[1]
data_list = read_lines_and_create_lists(filename)
# print(data_list)
new_list = process_list(data_list)
# print(new_list)

# Step 3: Opening a CSV file in write mode
out_name = filename.split('.vcf')[0] + '_counts_of_interest.csv'
with open(out_name, 'w', newline='') as file:
    # Step 4: Using csv.writer to write the list to the CSV file
    writer = csv.writer(file)
    writer.writerows(new_list) # Use writerows for nested list
