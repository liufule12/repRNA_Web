__author__ = 'aleeee'

import os
import shlex
import subprocess
import sys
import pickle
import math

from flask import request
from werkzeug.utils import secure_filename
import numpy as np
import matplotlib.pyplot as plt

import webserver.const as const


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in const.ALLOWED_EXTENSIONS


def create_user_fold(user_dir):
    if not os.path.exists(user_dir):
        cmd = 'mkdir ' + user_dir
        cmd_args = shlex.split(cmd)
        subprocess.Popen(cmd_args).wait()


def get_ind_names(filename):
    ind_names = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '>':
                ind_names.append(line.rstrip()[1:])

    return ind_names


def tran_args(form, mode):
    """transform the args to dict and add parameter k."""
    args = {}

    # Transform k.
    if mode == 'Triplet':
        args['k'] = 3
    elif mode == 'Kmer':
        args['k'] = int(form['k'])
    elif mode == 'PseSSC':
        args['k'] = int(form['n'])
    elif mode == 'PseDPC':
        args['k'] = int(form['d'])

    if mode in const.METHODS_PHYCHE_INDEX:
        args['k'] = 2

    # Transform lag.
    if mode in const.METHODS_LAG:
        args['lag'] = int(form['lag'])

    # Transform lamada, w.
    if mode in const.METHODS_LAMADA_W:
        args['lamada'] = int(form['lamada'])
        args['w'] = float(form['w'])

    return args


def save_file(filename, user_dir):
    try:
        rec_upload_file = request.files[filename]
        upload_file = secure_filename(rec_upload_file.filename)
        upload_file_path = user_dir + '/' + upload_file
        rec_upload_file.save(upload_file_path)
        return upload_file_path
    except:
        return None


def check_user_data(method, rec_data, form_args, input_file, write_file):
    if method == 'Kmer' or method == 'Triplet':
        check_res = write_form_file(rec_data, input_file, write_file, form_args['k'], 0, const.ALPHABET_RNA)
    elif method in const.METHODS_LAG:
        check_res = write_form_file(rec_data, input_file, write_file, form_args['k'], form_args['lag'], const.ALPHABET_RNA)
    elif method in const.METHODS_LAMADA_W:
        check_res = write_form_file(rec_data, input_file, write_file, form_args['k'], form_args['lamada'], const.ALPHABET_RNA)
    else:
        print("check_user_data RNA error!")

    return check_res


def is_alphabet(s, alphabet):
    """Check all element of line is under alphabet.
    """
    import re
    if alphabet == const.ALPHABET_RNA and re.search(r'[^ACGUacgu]', s) is not None:
        return False
    return True


def write_form_file(receive_data, filename, write_path, k, lamada, alphabet):
    """Receive the receive_data or file, then judge the data is FASTA form or not.
    Input: way_choice, the textarea file or upload file, write_path, LAMADA.

    Return (is_legal, count_seq, error_reason)

    TODO: optimize this ugly code.
    """
    if 0 == len(filename):
        lines = receive_data.split('\n')
    else:
        with open(filename) as f:
            lines = f.readlines()
    with open(write_path, 'w') as f:
        count_seq = 0
        len_seq = lamada  # Avoid the length of first seq is less than lamada.
        first_is_seq = True
        need_newline = False
        has_seq = True
        for e in lines:
            e = e.strip()
            if 0 == len(e):
                continue
            if '>' == e[0]:
                count_seq += 1
                # Sequence has not name.
                if len(e) <= 1:
                    return False, "Error, the sequence " + str(count_seq) + " has no name."
                # If the first_seq is not >, then error!
                if not first_is_seq:
                    return False, "Error, the sequence " + str(count_seq) + ": begin is not >."
                # If the sequence is not follow the >, then error!
                if not has_seq:
                    return False, "Error, the sequence " + str(count_seq - 1) + " has not sequence."
                # If the length of sequence is less than lamada, then error!
                if len_seq < lamada:
                    return False, "Error, the sequence " + str(count_seq - 1) + " parameter must be less than L-k, where L is the length of the query sequence and k is the length of the selected oligomer mode."
                if need_newline:
                    f.write('\n')
                need_newline = True
                f.write(e)
                f.write('\n')
                len_seq = 0
                first_is_seq = True
                has_seq = False
                continue
            e = e.upper()
            # If the sequence is not match the fasta form, then error!
            if 0 == count_seq and False == is_alphabet(e, alphabet):
                return False, "Error, the sequence " + str(count_seq + 1) + " has non " + alphabet + " characters."
            if not is_alphabet(e, alphabet):
                return False, "Error, the sequence " + str(count_seq) + " has non " + alphabet + " characters."
            if 0 == count_seq and True == first_is_seq:
                first_is_seq = False
            has_seq = True
            len_seq += len(e)
            f.write(e)

        # If the last sequence is not follow the >, then error!
        if not has_seq:
            return False, "Error, the sequence " + str(count_seq) + " has not sequence."
        # If there are only a sequence, but no >, then error!
        if count_seq < 1:
            return False, "Error, the sequence " + str(count_seq + 1) + " has not sequence name."
        # len_seq must be >= k.
        if len_seq < k:
            return False, "Error, the sequence " + str(count_seq) + " must be larger and equal to k, where k is the length of the selected oligomer mode."
        # If the length of sequence is less than lamada, then error!
        if len_seq < lamada:
            return False, "Error, the sequence " + str(count_seq) + " parameter must be less than L-k, where L is the length of the query sequence and k is the length of the selected oligomer mode."
        # PseSSC seq count must be less than 50.
        if (alphabet == const.ALPHABET_RNA) and count_seq > const.MAX_SEQ_NUM:
            return False, "Error, the sequence number must be less than " + str(const.MAX_SEQ_NUM) + "."

    return True, None, None, count_seq


def pse_process(method, args, input_file, ind_file, bracket_file, matched_file, vecs_file, user_dir):
    if method == 'Triplet':
        from pseALL.kmer import make_kmer_vector
        return make_kmer_vector(k=3, alphabet=const.ALPHABET_RNA, filename=input_file)
    elif method == 'Kmer':
        from pseALL.kmer import make_kmer_vector
        return make_kmer_vector(k=args['k'], alphabet=const.ALPHABET_RNA, filename=input_file)
    elif method in 'DAC':
        from pseALL.acc import acc
        return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
                   phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA)
    elif method in 'DCC':
        from pseALL.acc import acc
        return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
                   phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA, theta_type=2)
    elif method in 'DACC':
        from pseALL.acc import acc
        return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
                   phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA, theta_type=3)
    elif method == 'PC-PseDNC-General':
        from pseALL.pse import pseknc
        return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
                      phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA)
    elif method == 'SC-PseDNC-General':
        from pseALL.pse import pseknc
        return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
                      phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA, theta_type=2)
    elif method == 'PseSSC':
        generate_bracket_seq(input_file, bracket_file)
        psessc(args, bracket_file, vecs_file)
        mv_ps_file(bracket_file, user_dir)
        return read_tab_vecs(vecs_file)
    elif method == 'PseDPC':
        generate_bracket_seq(input_file, bracket_file)
        match_2st_has_name(bracket_file, matched_file)
        psedpc(args, matched_file, vecs_file)
        mv_ps_file(bracket_file, user_dir)
        return read_tab_vecs(vecs_file)
    else:
        print("pse_process error!")


def generate_bracket_seq(receive_file_path, bracket_file_path):
    """ This is a system command to generate bracket_seq file according receive_file. """
    cmd = "RNAfold <" + receive_file_path + " >" + bracket_file_path
    subprocess.Popen(cmd, shell=True).wait()


def psessc(args, bracket_file, vecs_file):
    """ Call java jar to generate PseSSC vecs."""
    cmd = "java -jar " + const.MODEL_PSESSC_PATH + " " + str(args['k']) + " " + str(args['lamada']) + " " + str(float(args['w'] * 10)) +  \
          " " + bracket_file + " " + vecs_file
    cmd_args = shlex.split(cmd)
    subprocess.Popen(cmd_args).wait()


def psedpc(args, bracket_file, vecs_file):
    """ Call java jar to generate PseDPC vecs."""
    cmd = "java -jar " + const.MODEL_PSEDPC_PATH + " " + str(args['k']) + " " + str(args['lamada']) + " " + str(float(args['w'] * 10)) +  \
          " " + bracket_file + " " + vecs_file
    cmd_args = shlex.split(cmd)
    subprocess.Popen(cmd_args).wait()


def match_2st_has_name(read_name, write_name):
    """
    Get matched sequence from 2st structure sequence and brackets.
    File per line is seq_name, seq_old and seq_matched included MFE.
    Input read file name, write file name.
    Return seq, seq is a list of triple(seq_name, seq_old, seq_bracket, seq_matched).
    """
    f_open = open(read_name)
    f_write = open(write_name, 'w')
    lines = f_open.readlines()
    len_lines = len(lines)
    seq = []

    for i in range(0, len_lines, 3):
        # Get seq_name.
        seq_name = ''
        for j in range(1, len(lines[i])):
            if ' ' != lines[i][j]:
                seq_name += lines[i][j]

        # Get the sequence_match sequence.
        sequence_old = lines[i + 1].strip()
        sequence_bracket = lines[i + 2].strip()
        stack_sequence = []
        stack_bracket = []
        dict_match = {}
        len_line = len(sequence_old)
        for j in range(0, len_line):

            if '.' == sequence_bracket[j]:
                dict_match[j] = '.'
            elif '(' == sequence_bracket[j]:
                stack_bracket.append('(')
                stack_sequence.append(j)
                dict_match[j] = '('
            elif ')' == sequence_bracket[j]:
                stack_bracket.pop()
                temp_order = stack_sequence.pop()
                dict_match[temp_order] = sequence_old[j]
                dict_match[j] = sequence_old[temp_order]
            else:
                print "match_2st error!!!!!!!!!!"
                return

        # Write seq_name, seq_old, seq_matched and include MFE.
        f_write.write(seq_name)
        f_write.write(lines[i + 1])
        str_match_values = ''.join(dict_match.values())
        f_write.write(str_match_values)
        # f_write.write(lines[i + 2][len_line:])
        f_write.write('\n')
        seq.append((seq_name.strip(), sequence_old, lines[i + 2][:len_line], str_match_values))
    return seq


def mv_ps_file(bracket_file, user_dir):
    """ Move the ps files (they are come from RNAfold command) to user fold."""
    with open(bracket_file) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i % 3 == 0:
                mv_file_name = line.rstrip()[1:] + "_ss.ps"
                cmd = "mv " + os.getcwd() + "/" + mv_file_name + " " + user_dir + "/" + mv_file_name
                subprocess.Popen(cmd, shell=True).wait()


def read_tab_vecs(read_file):
    with open(read_file) as f:
        lines = f.readlines()
        vecs = []

        for line in lines:
            vecs.append(line.rstrip().split())

    return vecs


def write_tab(mode, args, _vecs, vecs_name, write_file):
    """Write the vectors into disk in tab format."""
    with open(write_file, 'w') as f:
        # Write the parameters.
        f.write("Data type: RNA sequences\n")
        f.write("Mode: " + mode + '\n')
        if mode == 'PseSSC':
            f.write('N: ' + str(args['k']) + '\n')
        elif mode == 'PseDPC':
            f.write('D: ' + str(args['k']) + '\n')
        else:
            f.write('K: ' + str(args['k']) + '\n')
        if args['props']:
            f.write('Properties: ' + str(args['props']) + '\n')
        if args['ext_ind']:
            f.write("User-defined properties: " + str(args['ext_ind']) + '\n')
        if 'lamada' in args:
            f.write('Lambda: ' + str(args['lamada']) + '\n')
        elif 'lag' in args:
            f.write('Lag: ' + str(args['lag']) + '\n')
        if 'w' in args:
            f.write('W: ' + str(args['w']) + '\n')
        f.write('\n')

        # Write the results.
        for ind, vec in enumerate(_vecs):
            f.write('>' + str(vecs_name[ind]))
            f.write('\n')
            f.write(str(vec[0]))
            for val in vec[1:]:
                f.write('\t' + str(val))
            f.write('\n')


def heatmap(data, write_file, args):
    len_data = len(data)
    line_num = int(math.ceil(math.sqrt(len_data)))

    if args['mode'] in const.METHODS_LAMADA_W:
        # Process the data.
        tar_data = [0] * 2
        spl_num = len_data - args['lamada']
        per_num = spl_num / line_num
        tar_data[0] = get_fir_data(data[:spl_num], per_num, line_num)
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
        # Process the data.
        per_num = len_data / line_num
        tar_data = get_fir_data(data, per_num, line_num)

        # Plot the pic.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if args['mode'] in const.METHODS_ACC:
            ax.set_title("Distribution of sequence order correlation values")
        else:
            ax.set_title("Distribution of sequence composition values")

        ax.axes.get_yaxis().set_visible(False)
        cax = ax.imshow(tar_data, interpolation='nearest')
        cbar = fig.colorbar(cax)
        fig.savefig(write_file)


def get_fir_data(data, per_num, line_num):
    """Change the vector data to matrix data for ploting."""
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

        # Supply the none data in last line.
        fir_data = temp_data
        try:
            if len(fir_data[-1]) != len(fir_data[-2]):
                temp_min = min(data)
                for i in range(len(fir_data[-1]), len(fir_data[-2])):
                    fir_data[-1].append(temp_min)
        except IndexError:
            pass

    return fir_data