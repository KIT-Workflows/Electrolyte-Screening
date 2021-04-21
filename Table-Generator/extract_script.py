import numpy as np
import os, sys, re, csv, fnmatch, yaml

def separate_string_number(string):
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i

        previous_character = i

        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups

def find_parameter(file, Search_parameter):
    parameter_var = []
    nums = re.compile(r"[+-]?\d+(?:\.\d+)?")  # Regular expression
    with open(file, "r") as file1:
        for line in file1.readlines():
            if Search_parameter in line:
                idx1 = line.find('"')
                idx2 = line.find('"', idx1 + 1)
                field = line[idx1 + 1:idx2]
                array = nums.findall(field)
                # print(array)
                # print(Search_parameter)
                
                if len(array) > 0:
                    if Search_parameter == 'ENCUT':
                        parameter_var = np.around(float(array[0]), 3)
                    else:    
                        parameter_var = np.around(float(array[len(array)-1]), 3)
                else:
                    continue
        return parameter_var


def find_var1(file, Search_string):
    with open(file, "r") as file1:
        for line in file1.readlines():
            if Search_string in line:
                idx1 = line.find('"')
                idx2 = line.find('"', idx1 + 2)
                field = line[idx1 + 1:idx2 - 1]
        return str.strip(field.rpartition(' = ')[2])

def Filter(string, substr): 
	return[str1 for str1 in string if
			any(sub in str1 for sub in substr)] 


def round_float_list(float_list, decimal_points):
    float_list = [round(float(item),decimal_points) for item in float_list]
    return float_list

if __name__ == '__main__':
    # Input data
    decimal_points = 4
    with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)
    
    Search_parameter_list = []
    for ii in range(len(wano_file["Search-Parameters"])):
        Search_parameter_list.append(wano_file["Search-Parameters"][ii]["var"])

    Search_file = wano_file["Search-in-file"]

    Search_file_list = []
    Table_var = []

    for file1 in os.listdir():
        if file1.endswith(Search_file):
            Search_file_list.append(file1)

    temp_len = 0
    for file in Search_file_list:
        for parameter_var in Search_parameter_list:
            with open(file) as f_file:
                n_file = yaml.full_load(f_file)
            temp_var1 = n_file[parameter_var]
            Table_var.append(temp_var1)
        # else:
        #     for parameter_var in Search_parameter_list:
        #         temp_var1 = find_parameter(file, parameter_var)
        #         Table_var.append(temp_var1)
    
    Table_var = np.array(Table_var)
    Table_var = Table_var.reshape(len(Search_file_list), len(Search_parameter_list))
    Sorted_Table = Table_var[np.argsort(Table_var[:, 1])]

    Table_dict = {}
    for ii in range(len(Search_parameter_list)):
        Table_dict[Search_parameter_list[ii]] =  round_float_list(Sorted_Table[:,ii],decimal_points) 

    try:   
        with open("Table-var-dict.yml",'w') as out:
            yaml.dump(Table_dict, out,default_flow_style=False)
    
        with open('Table-var', 'w') as f: 
            # using csv.writer method from CSV package 
            write = csv.writer(f) 
            # Columns label
            write.writerow(Search_parameter_list) 
            # Columns values
            write.writerows(Sorted_Table)
    except IOError:
        print("I/O error") 

    if wano_file["Delete-Files"]:
        for file in Search_file_list:
            os.remove(file)