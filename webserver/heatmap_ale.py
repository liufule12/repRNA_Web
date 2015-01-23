__author__ = 'aleeee'

import sys
import pickle
import math
import numpy as np
import matplotlib.pyplot as plt
import webserver.const as const


def get_fir_data(data, per_num):
    fir_data = []
    print(len(data))
    if per_num == 0:
        fir_data = [data]
    elif per_num > 0:
        temp_data = []
        for ind in range(per_num + 1):
            temp_line = data[ind*line_num:(ind+1)*line_num]
            if temp_line:
                temp_data.append(temp_line)
        fir_data = temp_data
        try:
            if len(fir_data[-1]) != len(fir_data[-2]):
                temp_min = min(data)
                for i in range(len(fir_data[-1]), len(fir_data[-2])):
                    fir_data[-1].append(temp_min)
        except IndexError:
            pass

    return fir_data


if __name__ == '__main__':
    print(sys.argv)

    read_file = sys.argv[1]
    write_file = sys.argv[2]
    alp_num = const.RNA_ALP_NUM
    method = sys.argv[3]
    k = int(sys.argv[4])
    try:
        lamada = int(sys.argv[5])
    except IndexError:
        pass

    with open(read_file, 'rb') as vis_f:
        data = pickle.load(vis_f)
        len_data = len(data)

        if method in const.METHODS_PSE:
            # Process the data.
            line_num = 64
            tar_data = [0] * 2
            spl_num = alp_num ** k
            per_num = spl_num / line_num
            tar_data[0] = get_fir_data(data[:spl_num], per_num)
            tar_data[1] = [data[spl_num:]]

            # Plot the pic.
            fig, axes = plt.subplots(nrows=2, ncols=1)

            ax = axes.flat[0]
            im = ax.imshow(tar_data[0], interpolation='nearest')
            ax.set_title("Distribution of sequence composition values")
            ax.axes.get_yaxis().set_visible(False)

            bx = axes.flat[1]
            im = bx.imshow(tar_data[1], interpolation='nearest')
            bx.set_title("Distribution of sequence order correlation values")
            bx.axes.get_yaxis().set_visible(False)

            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(im, cax=cbar_ax)

            fig.savefig(write_file)
        else:
            line_num = 64
            per_num = len_data / line_num
            tar_data = get_fir_data(data, per_num)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            if method in const.METHODS_ACC:
                ax.set_title("Distribution of sequence order correlation values")
            else:
                ax.set_title("Distribution of sequence composition values")

            ax.axes.get_yaxis().set_visible(False)
            cax = ax.imshow(tar_data, interpolation='nearest')
            cbar = fig.colorbar(cax)
            fig.savefig(write_file)