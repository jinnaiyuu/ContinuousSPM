import sys
import numpy as np
import csv

def binarize_median(input_filename, output_filename, output_header):
    data = np.loadtxt(input_filename, delimiter=",")
    col_names = []
    col_thrs = []
    for j in range(data.shape[1]):
        col_names.append("[<="+str(np.median(data[:,j]))+"]")
        col_names.append("[>"+str(np.median(data[:,j]))+"]")
        col_thrs.append(np.median(data[:,j]))
    data_trans = []
    for i_x in range(data.shape[0]):
        x = []
        idx = 0
        for j in range(data.shape[1]):
            if data[i_x][j] <= col_thrs[j]:
                x.append(idx)
            idx += 1
            if data[i_x][j] > col_thrs[j]:
                x.append(idx)
            idx += 1
        data_trans.append(x)
    f = open(output_filename, 'w')
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_trans)
    f.close()
    f = open(output_header, 'w')
    writer = csv.writer(f, delimiter=' ')
    writer.writerow(col_names)
    f.close()

def binarize_interordinal(input_filename, output_filename, output_header, digit):
    data = np.loadtxt(input_filename, delimiter=",")
    data = np.around(data, decimals=digit)
    col_names = []
    for j in range(data.shape[1]):
        for i in range(data.shape[0]):
            col_names.append("[<="+str(data[i][j])+"]")
            col_names.append("[>="+str(data[i][j])+"]")
    data_trans = []
    for i_x in range(data.shape[0]):
        x = []
        idx = 0
        for j in range(data.shape[1]):
            for i in range(data.shape[0]):
                if data[i_x][j] <= data[i][j]:
                    x.append(idx)
                idx += 1
                if data[i_x][j] >= data[i][j]:
                    x.append(idx)
                idx += 1
        data_trans.append(x)
    f = open(output_filename, 'w')
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_trans)
    f.close()
    f = open(output_header, 'w')
    writer = csv.writer(f, delimiter=' ')
    writer.writerow(col_names)
    f.close()

def binarize_interval(input_filename, output_filename, output_header, digit):
    data = np.loadtxt(input_filename, delimiter=",")
    data = np.around(data, decimals=digit)
    col_names = []
    for j in range(data.shape[1]):
        for i in range(data.shape[0] - 1):
            for i2 in range(i + 1, data.shape[0]):
                col_names.append("("+str(min(data[i][j], data[i2][j]))+","+str(max(data[i][j], data[i2][j]))+"]")
    data_trans = []
    for i_x in range(data.shape[0]):
        x = []
        idx = 0
        for j in range(data.shape[1]):
            for i in range(data.shape[0] - 1):
                for i2 in range(i + 1, data.shape[0]):
                    if min(data[i][j], data[i2][j]) < data[i_x][j] <= max(data[i][j], data[i2][j]):
                        x.append(idx)
                    idx += 1
        data_trans.append(x)
    f = open(output_filename, 'w')
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_trans)
    f.close()
    f = open(output_header, 'w')
    writer = csv.writer(f, delimiter=' ')
    writer.writerow(col_names)
    f.close()


digit = 5
## get arguments
argv = sys.argv
argc = len(argv)
## check the nubmer of args
if argc < 5:
    print("Error: lack of args")
    print("Usage: python " + argv[0] + " <input_file_name> <output_file_name> <output_header_name> <binarization_type>")
    quit()

## get arguments
input_file = argv[1]
output_file = argv[2]
output_header = argv[3]
binarization_type = int(argv[4])
if argc == 6:
    digit = int(argv[5])

## perform binarization
if binarization_type == 1:
    binarize_median(input_file, output_file, output_header)
elif binarization_type == 2:
    binarize_interordinal(input_file, output_file, output_header, digit)
elif binarization_type == 3:
    binarize_interval(input_file, output_file, output_header, digit)
else:
    binarize_median(input_file, output_file, output_header)
