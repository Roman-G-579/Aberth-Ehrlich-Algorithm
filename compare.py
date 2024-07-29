import re
#######################################
#If you are using this the roman
# the result.txt file contain python output [(complex number) , (complex number) , (complex number]
# the result2.txt file contain scheme output (complex_number complex_number complex_number)
# in this exact format 
######################################
# Function to extract complex numbers from result.txt
def extract_complex_numbers_result(txt):
    complex_pattern = re.compile(r'\(([^)]+)\)')
    complex_numbers = []
    for match in complex_pattern.findall(txt):
        complex_numbers.append(complex(match))
    return complex_numbers

# Function to extract complex numbers from result2.txt
def extract_complex_numbers_result2(txt):
    complex_pattern = re.compile(r'\(([^)]+\))')
    complex_numbers = []
    for match in complex_pattern.findall(txt.replace('i', 'j'))[0].replace(')', '').split(' '):
        print(match)
        complex_numbers.append(complex(match))
    return complex_numbers

# Read the contents of the files
with open('result.txt', 'r', encoding='utf-16') as file1:
    content_result = file1.read()

with open('result2.txt', 'r', encoding='utf-16') as file2:
    content_result2 = file2.read()

# Extract complex numbers from both files
complex_numbers_result = extract_complex_numbers_result(content_result)
complex_numbers_result2 = extract_complex_numbers_result2(content_result2)

# Find common complex numbers
common_complex_numbers = set(complex_numbers_result).intersection(complex_numbers_result2)

# find numbers in result but not in result2
numbers_in_result_not_in_result2 = set(complex_numbers_result).difference(complex_numbers_result2)

# find numbers in result2 but not in result
numbers_in_result2_not_in_result = set(complex_numbers_result2).difference(complex_numbers_result)

# Display the results
print("Numbers in result.txt but not in result2.txt:", numbers_in_result_not_in_result2)
print("Numbers in result2.txt but not in result.txt:", numbers_in_result2_not_in_result)


