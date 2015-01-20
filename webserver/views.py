from webserver import app

from flask import Flask, render_template, redirect, request

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

        # Check and transform the form args.
        check_res = _method.tran_args(request.form, mode)
        if check_res[0] is False:
            return render_template('result.html', er_info=(True, check_res[1]))
        else:
            form_args = check_res[1]
        print("Args is ok.")

        return "Process in main completed."

