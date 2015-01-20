__author__ = 'aleeee'

import os
import shlex
import subprocess

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
            if line == '>':
                ind_names.append(line.rstrip()[1:])

    return ind_names


def tran_args(form, mode):
    """Transform the args type."""
    print("Tran_args", form, mode)
    args = {}

    # Transform k.
    if mode == 'Kmer':
        k = int(form['k'])
        if k < 1 or k > 6:
            return False, "Error: The value of k must be an integer and larger than 0 and small than 7."
        else:
            args['k'] = k
    if mode in const.METHODS_PHYCHE_INDEX:
        args['k'] = 2

    # Transform lag.
    if mode in const.METHODS_LAG:
        try:
            lag = int(form['lag'])
        except ValueError:
            return False, "Error: The value of lag must be an integer larger than 0 and smaller than " \
                          "the length of the longest input sequence - lag."

        if lag < 0:
            return False, "Error: The value of lag must be an integer larger than 0 and smaller than " \
                          "the length of the longest input sequence - lag."
        else:
            args['lag'] = lag

    # Transform n.
    if mode == 'PseSSC':
        try:
            n = int(form['n'])
        except ValueError:
            return False, "Error: n error information."

        if n < 0:
            return False, "Error: The value of lag must be an integer larger than 0 and n error information."
        else:
            args['n'] = n

    # Transform d.
    if mode == 'PseDPC':
        try:
            d = int(form['d'])
        except ValueError:
            return False, "Error: d error information."

        if d < 0:
            return False, "Error: The value of lag must be an integer larger than 0 and n error information."
        else:
            args['d'] = d

    # Transform lamada, w.
    if mode in const.METHODS_LAMADA_W:
        try:
            lamada = int(form['lamada'])
        except ValueError:
            return False, "Error: The value of lambda must be an integer larger than 0 and smaller than " \
                          "the length of the longest input sequence - lambda."

        try:
            w = float(form['w'])
        except ValueError:
            return False, "Error: The value of w must be no less than 0 and no larger than 1."

        if lamada < 0:
            return False, "Error: The value of lambda must be an integer larger than 0 and smaller than " \
                          "the length of the longest input sequence - lambda."
        elif w < 0 or w > 1:
            return False, "Error: The value of w must be no less than 0 and no larger than 1."
        else:
            args['lamada'] = lamada
            args['w'] = w

    return True, args


def check_user_data(category, method, request, form_args, input_file, write_file):
    print("form_args:", form_args)
    if category == const.DNA:
        if method == 'Kmer' or method == 'RevKmer':
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, 0, const.ALPHABET_DNA)
        elif method in const.DNA_PSE:
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, form_args['lamada'], const.ALPHABET_DNA)
        elif method in const.DNA_ACC:
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, form_args['lag'], const.ALPHABET_DNA)
        else:
            print("check_user_data DNA error!")
    elif category == const.RNA:
        if method == 'Kmer':
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, 0, const.ALPHABET_RNA)
        elif method in const.RNA_PSE:
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, form_args['lamada'], const.ALPHABET_RNA)
        elif method in const.RNA_ACC:
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, form_args['lag'], const.ALPHABET_RNA)
        else:
            print("check_user_data RNA error!")
    elif category == const.PROTEIN:
        if method == 'Kmer':
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, 0, const.ALPHABET_PROTEIN)
        elif method in const.PROTEIN_PSE:
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, form_args['lamada'], const.ALPHABET_PROTEIN)
        elif method in const.PROTEIN_ACC:
            check_res = write_form_file(request.form['rec_data'], input_file, write_file, form_args['lag'], const.ALPHABET_PROTEIN)
        else:
            print("check_user_data PROTEIN error!")
    else:
        print("check_user_data error")

    return check_res


def is_alphabet(s, alphabet):
    """Check all element of line is under alphabet.
    """
    import re

    if alphabet == const.ALPHABET_DNA and re.search(r'[^ACGTacgt]', s) is not None:
        return False
    elif alphabet == const.ALPHABET_RNA and re.search(r'[^ACGUacgu]', s) is not None:
        return False
    elif alphabet == const.ALPHABET_PROTEIN and re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', s) is not None:
        return False
    return True


def write_form_file(receive_data, filename, write_path, lamada, alphabet):
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
        # If the length of sequence is less than lamada, then error!
        if len_seq < lamada:
            return False, "Error, the sequence " + str(count_seq) + " parameter must be less than L-k, where L is the length of the query sequence and k is the length of the selected oligomer mode."
        # PseSSC seq count must be less than 50.
        if (alphabet == const.ALPHABET_DNA or const.ALPHABET_RNA) and count_seq > const.MAX_DNA_NUM:
            return False, "Error, the sequence number must be less than " + str(const.MAX_DNA_NUM) + "."
        # Seqs numbers must be less than 10.
        if alphabet == const.ALPHABET_PROTEIN and count_seq > const.MAX_PROTEIN_NUM:
            return False, "Error, the sequence number must be less than " + str(const.MAX_PROTEIN_NUM) + "."
    return True, None, None, count_seq


# def pse_process(category, method, args, input_file, ind_file):
#     print("Pse_Process args", category, method, args, input_file)
#
#     if category == const.DNA:
#         if method == 'Kmer':
#             from pseALL.kmer import make_kmer_vector
#             return make_kmer_vector(k=args['k'], alphabet=const.ALPHABET_DNA, filename=input_file)
#         elif method == 'RevKmer':
#             from pseALL.kmer import make_kmer_vector
#             return make_kmer_vector(k=args['k'], alphabet=const.ALPHABET_DNA, filename=input_file, revcomp=True)
#         elif method == 'PseDNC' or \
#                 method == 'PC-PseDNC-General' or method == 'PC-PseTNC-General':
#             from pseALL.pse import pseknc
#             return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
#                           phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_DNA)
#         elif method == 'SC-PseDNC-General' or method == 'SC-PseTNC-General':
#             from pseALL.pse import pseknc
#             return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
#                           phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_DNA, theta_type=2)
#         elif method == 'PseKNC':
#             from pseALL.pse import ipseknc
#             return ipseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
#                            phyche_list=args['props'], alphabet=const.ALPHABET_DNA)
#         elif method == 'DAC' or method == 'TAC':
#             from pseALL.acc import acc
#             return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                        phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_DNA)
#         elif method == 'DCC' or method == 'TCC':
#             from pseALL.acc import acc
#             return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                        phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_DNA, theta_type=2)
#         elif method == 'DACC' or method == 'TACC':
#             from pseALL.acc import acc
#             return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                        phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_DNA, theta_type=3)
#     elif category == const.RNA:
#         if method == 'Kmer':
#             from pseALL.kmer import make_kmer_vector
#             return make_kmer_vector(k=args['k'], alphabet=const.ALPHABET_RNA, filename=input_file)
#         elif method == 'PC-PseDNC-General':
#             from pseALL.pse import pseknc
#             return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
#                           phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA)
#         elif method == 'SC-PseDNC-General':
#             from pseALL.pse import pseknc
#             return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
#                           phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA, theta_type=2)
#         elif method in const.RNA_ACC:
#             from pseALL.acc import acc
#             return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                        phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_RNA)
#     elif category == const.PROTEIN:
#         if method == 'Kmer':
#             from pseALL.kmer import make_kmer_vector
#             return make_kmer_vector(k=args['k'], alphabet=const.ALPHABET_PROTEIN, filename=input_file)
#         elif method in const.PROTEIN_PSE:
#             from pseALL.pse import pseknc
#             return pseknc(input_data=open(input_file), k=args['k'], w=args['w'], lamada=args['lamada'],
#                           phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_PROTEIN)
#         elif method in const.PROTEIN_ACC:
#             from pseALL.acc import acc
#             if method == 'AC':
#                 return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                            phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_PROTEIN)
#             elif method == 'CC':
#                 return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                            phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_PROTEIN, theta_type=2)
#             elif method == 'ACC':
#                 return acc(input_data=open(input_file), k=args['k'], lag=args['lag'],
#                            phyche_list=args['props'], extra_index_file=ind_file, alphabet=const.ALPHABET_PROTEIN, theta_type=3)


def write_tab(category, method, args, _vecs, vecs_name, write_file):
    """Write the vectors into disk in tab format."""
    with open(write_file, 'w') as f:
        # Write the parameters.
        if category == const.PROTEIN:
            f.write("Data type: " + 'protein sequences\n')
        elif category == const.DNA or category == const.RNA:
            f.write("Data type: " + category + ' sequences\n')
        f.write("Method: " + method + '\n')
        if 'k' in args:
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