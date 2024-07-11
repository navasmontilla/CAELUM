import os
import math 
import argparse                   
import numpy as np             
import matplotlib.pyplot as plt 
from utils import modify_header_file, write_config, write_initial, write_equilibrium, backup_file, restore_file, backup_fileL2, restore_fileL2, compile_program, run_program, initialize_variables, read_data_euler, singleRP, caseBubble, caseLinear, ordersLinear

# Set up argument parsing
parser = argparse.ArgumentParser(description="Run test cases and print results.")
parser.add_argument('arg1', type=int, nargs='?', default=4, help="N threads")
parser.add_argument('arg2', type=int, nargs='?', default=0, help="Reconstruction method")
parser.add_argument('arg3', type=int, nargs='?', default=5, help="Order")
args = parser.parse_args()

script_dir = os.path.dirname(os.path.abspath(__file__))
folder_lib = os.path.join(script_dir, "../lib")
backup_fileL2(folder_lib+'/definitions.h')
modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', args.arg1)
modify_header_file(folder_lib+'/definitions.h', 'TYPE_REC', args.arg2)
ord=args.arg3

# Define colors for the output
bold_green = "\033[1;32m"
bold_red = "\033[1;31m"
bold_blue = "\033[1;34m"
reset = "\033[0m"

# Define the tests
test_names = ["Convergence test (linear)","RP1 (equilibrium)", "RP2 (sod-shock)", "RP3", "RP4 (2-component)", "Bubble collision"]

# Run the Linear transport test
l = ordersLinear()

# Run RP tests
rp = np.zeros(4)
for i in range(1, 5):
    a = singleRP(i,ord)
    rp[i-1] = a

# Run the Bubble collision case test
b = caseBubble(ord)

restore_fileL2(folder_lib+'/definitions.h')

# Store test results
results = []

results.append((test_names[0], l))

for i in range(1, 5):
    results.append((test_names[i], rp[i-1]))

# Include the Bubble test result
results.append((test_names[5], b))

# Print the results in a formatted table
print(f"{bold_blue}{'Test Name':<30}{'Result':<10}{reset}")
print(f"{bold_blue}{'-'*40}{reset}")

# Count the number of passed and failed tests
passed_count = 0
failed_count = 0

for test_name, result in results:
    if result == 1:
        print(f"{test_name:<30}{bold_green}{'Passed':<10}{reset}")
        passed_count += 1
    else:
        print(f"{test_name:<30}{bold_red}{'Failed':<10}{reset}")
        failed_count += 1

# Print a summary of the results
print(f"\n{bold_blue}{'Summary':<30}{'Count':<10}{reset}")
print(f"{bold_blue}{'-'*40}{reset}")
print(f"{'Passed tests:':<30}{bold_green}{passed_count:<10}{reset}")
print(f"{'Failed tests:':<30}{bold_red}{failed_count:<10}{reset}")
print(f"{'-'*40}")
print(f"{'Total tests:':<30}{passed_count + failed_count:<10}")