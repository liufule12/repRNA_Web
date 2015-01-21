import os
import time

from webserver import app
from flask import Flask, render_template, redirect, request

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

        # Save the upload file.
        rec_upload_file = _method.save_file('upload_data', user_dir)
        rec_ind_file = _method.save_file('upload_ind', user_dir)
        print(rec_upload_file, rec_ind_file)
        print("The user upload file is ok.")

        return "Process in main completed."


@app.route("/test/")
def test():
    return render_template("test.html")