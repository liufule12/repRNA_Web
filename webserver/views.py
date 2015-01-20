import os
import time

from webserver import app
from flask import Flask, render_template, redirect, request
from werkzeug.utils import secure_filename
from werkzeug.exceptions import RequestEntityTooLarge

import webserver._method as _method


@app.route('/')
@app.route('/home/')
def home():
    return render_template("home.html")


@app.route('/server/')
def server():
    return redirect("/RNA/kmer/")


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
        print("Args is ok.", form_args)

        # Create the user fold.
        user_ip_time = request.remote_addr + '_' + str(time.time())
        user_dir = os.getcwd() + '/webserver/static/temp/' + user_ip_time
        _method.create_user_fold(user_dir)
        print("The user fold is ok.")

        # Deal with the upload data file.
        rec_upload_file = request.files['upload_data']
        print(rec_upload_file)
        if rec_upload_file and _method.allowed_file(rec_upload_file.filename):
            upload_file = secure_filename(rec_upload_file.filename)
            upload_file_path = user_dir + '/' + upload_file
            rec_upload_file.save(upload_file_path, buffer_size=1)
        elif rec_upload_file and not _method.allowed_file(rec_upload_file.filename):
            return render_template("result.html",
                                   er_info=(True, "Sorry, the upload file must be txt or fasta file."))
        else:
            upload_file = None
        print("The user upload data file is ok.")

        return "Process in main completed."


@app.route("/test/")
def test():
    return render_template("test.html")