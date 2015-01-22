__author__ = 'aleeee'

import os
import shlex
import subprocess

import webserver.const as const
from flask import request
from werkzeug.utils import secure_filename


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
            if line == '>':
                ind_names.append(line.rstrip()[1:])

    return ind_names


def tran_args(form, mode):
    """transform the args to dict and add parameter k."""
    print("Tran_args", form, mode)
    args = {}

    # Transform k.
    if mode == 'Kmer':
        args['k'] = int(form['k'])
    if mode in const.METHODS_PHYCHE_INDEX:
        args['k'] = 2

    # Transform lag.
    if mode in const.METHODS_LAG:
        args['lag'] = int(form['lag'])

    # Transform n.
    if mode == 'PseSSC':
        args['n'] = int(form['n'])

    # Transform d.
    if mode == 'PseDPC':
        args['d'] = int(form['d'])

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
    print("form_args:", method, form_args, form_args['k'])
    if method == 'Kmer':
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
            # print e
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


def pse_process(method, args, input_file, ind_file):
    print("Pse_Process args", method, args, input_file)

    if method == 'Kmer':
        from pseALL.kmer import make_kmer_vector
        return make_kmer_vector(k=args['k'], alphabet=const.ALPHABET_RNA, filename=input_file)
    elif method in const.METHODS_LAG:
        from pseALL.acc import acc
        return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
                   phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA)
    elif method == 'PC-PseDNC-General':
        from pseALL.pse import pseknc
        return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
                      phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA)
    elif method == 'SC-PseDNC-General':
        from pseALL.pse import pseknc
        return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
                      phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA, theta_type=2)
    elif method == 'PseSSC':
        pass
    elif method == 'PseDPC':
        pass
    else:
        print("pse_process error!")


def write_tab(mode, args, _vecs, vecs_name, write_file):
    """Write the vectors into disk in tab format."""
    with open(write_file, 'w') as f:
        # Write the parameters.
        f.write("Data type: RNA sequences\n")
        f.write("Mode: " + mode + '\n')
        if 'k' in args:
            f.write('K: ' + str(args['k']) + '\n')
        if 'n' in args:
            f.write('N: ' + str(args['k']) + '\n')
        if 'd' in args:
            f.write('D: ' + str(args['k']) + '\n')
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