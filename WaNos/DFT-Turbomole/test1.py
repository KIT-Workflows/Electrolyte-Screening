def humo_lumo_number_orbt(file_name, string_to_search):
    #Searching for the given string in file along with its line numbers
    line_number = 0
    list_of_results = []
    # Opening the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Reading all lines in the file one by one by iterating the file
        for line in read_obj:
            # checking each line, if the line contains the string
            if string_to_search in line:
                # If it contains the string, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))
            line_number += 1

    # Return list of tuples containing line numbers and lines where string is found

    with open(file_name, 'r') as file:
        content = file.readlines()

    line_number = 1 + list_of_results[0][0]

    homo_line = int(content[line_number].split()[2])
    lumo_line = int(content[line_number+1].split()[2])
    
    return homo_line, lumo_line 
l1, l2 = humo_lumo_number_orbt('eiger.out', 'HOMO-LUMO Separation')
print(l1,l2)
