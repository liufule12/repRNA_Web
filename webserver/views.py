import os
import shlex
import subprocess
import time
import pickle

from webserver import app
from flask import Flask, render_template, redirect, request
from werkzeug.exceptions import RequestEntityTooLarge

import webserver._method as _method
import webserver.pseALL.util as util
import webserver.const as const


@app.route('/')
@app.route('/home/')
def home():
    return render_template("home.html")


@app.route('/server/')
def server():
    return redirect("/RNA/Kmer/")


@app.route('/tutorial.html')
def tutorial():
    return render_template("tutorial.html")


@app.route('/doc/')
def doc():
    return render_template("doc.html")


@app.route('/download/')
def download():
    return render_template("download.html")


@app.route('/citation/')
def citation():
    return render_template("citation.html")


@app.route('/contact/')
def contact():
    return render_template("contact.html")


@app.route('/RNA/<mode>/', methods=['GET', 'POST'])
def main(mode):
    if request.method == 'GET':
        return render_template("RNA.html", mode=mode)
    if request.method == 'POST':
        print("request.form", request.form)
        print("request.files", request.files)

        # Transform the form args and add parameter k.
        form_args = _method.tran_args(request.form, mode)
        form_args['props'] = [e for e in request.form if e not in const.ARGS]
        print("Args is ok.", form_args)

        # Create the user fold.
        user_ip_time = request.remote_addr + '_' + str(time.time())
        user_dir = os.getcwd() + '/webserver/static/temp/' + user_ip_time
        _method.create_user_fold(user_dir)
        print("The user fold is ok.")

        # Save the upload file.
        data_file_path = _method.save_file('upload_data', user_dir)
        ind_file_path = _method.save_file('upload_ind', user_dir)
        form_args['ext_ind'] = []
        if ind_file_path is not None:
            form_args['ext_ind'] = _method.get_ind_names(ind_file_path)
        print(form_args['ext_ind'])
        print("The user upload file is ok.")

        # Check the user data and write the data file into user directory.
        rec_data = request.form['rec_data']
        if 0 != len(rec_data):
            data_file = ""
        else:
            data_file = data_file_path

        input_file = user_dir + '/' + 'input.txt'
        check_res = _method.check_user_data(method=mode, rec_data=rec_data, form_args=form_args,
                                            input_file=data_file, write_file=input_file)
        if check_res[0] is False:
            return render_template("result.html", er_info=(True, check_res[1]))
        print("rec_data is ok.")

        # Get sequences names.
        seqs = util.get_data(input_data=open(input_file), alphabet="RNA", desc=True)
        names = [seq.name for seq in seqs]
        print(names)
        print("seq names is ok.")

        # Process.
        try:
            bracket_file = user_dir + '/' + 'bracket.txt'
            matched_file = user_dir + '/' + 'matched.txt'
            vecs_file = user_dir + '/' + 'vecs.txt'
            res = _method.pse_process(method=mode, args=form_args, input_file=input_file, ind_file=ind_file_path,
                                      bracket_file=bracket_file, matched_file=matched_file, vecs_file=vecs_file)
        except:
            if ind_file_path is not None:
                return render_template("result.html",
                                       er_info=(True, "The physicochemical indices file format error."))
            raise

        # Write the res in TAB format.
        write_file = user_dir + '/res.txt'
        download_path = 'temp/' + user_ip_time + '/res.txt'
        _method.write_tab(mode=mode, args=form_args, vecs_name=names, _vecs=res, write_file=write_file)

        print(res)
        return render_template('result.html', er_info=(False, None), res=res, mode=mode, args=form_args, names=names,
                               write_file=download_path, user_ip_time=user_ip_time)


@app.route("/visual/<user_ip_time>/<ind>/")
def visual(user_ip_time, ind):
    ind = int(ind)
    args = {}
    lines = []
    user_dir = os.getcwd() + '/webserver/static/temp/' + user_ip_time
    with open(user_dir + '/res.txt') as f:
        lines = f.readlines()
    for line_ind, line in enumerate(lines):
        if line != '\n':
            line = line.split(': ')
            if line[0] == 'Data type':
                args['category'] = line[1].rstrip().split(' ')[0]
            elif line[0] == 'Mode':
                args['mode'] = line[1].rstrip()
            elif line[0] == 'K':
                args['k'] = int(line[1])
            elif line[0] == 'Lambda':
                args['lamada'] = int(line[1])
            continue
        lines = lines[line_ind+1:]
        break

    # Find the visualization line.
    vec_name = lines[ind * 2].rstrip()
    vis_line = lines[ind * 2 + 1].rstrip().split('\t')
    vis_line = [float(e) for e in vis_line]

    # Plot heatmap.
    vis_path = user_dir + '/vis.txt'
    tar_vis_path = user_dir + '/vis.jpg'
    with open(vis_path, 'wb') as f:
        pickle.dump(vis_line, f)
    try:
        cmd = 'python webserver/heatmap_ale.py ' + vis_path + ' ' + tar_vis_path + ' ' +\
              args['mode'] + ' ' + str(args['k']) + ' ' + str(args['lamada'])
    except KeyError:
        cmd = 'python webserver/heatmap_ale.py ' + vis_path + ' ' + tar_vis_path + ' ' +\
              args['mode'] + ' ' + str(args['k'])
    cmd_args = shlex.split(cmd)
    subprocess.Popen(cmd_args).wait()

    jpg_path = 'temp/' + user_ip_time + '/vis.jpg'
    return render_template("visualization.html", jpg_path=jpg_path, ind=ind, vec=vis_line, vec_name=vec_name, args=args)


@app.route("/test/")
def test():
    return render_template("test.html")