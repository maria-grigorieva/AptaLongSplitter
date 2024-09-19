from datetime import datetime

import pandas as pd
from flask import Flask, render_template, request, url_for, redirect, session, jsonify
import requests
from flask_bootstrap import Bootstrap5
from flask_wtf import FlaskForm, CSRFProtect
from plotly.subplots import make_subplots
from werkzeug.utils import secure_filename
from flask_wtf.file import FileField, FileRequired
from wtforms.validators import DataRequired, Length
from wtforms import StringField, FloatField, SubmitField, FieldList, FormField, Form, SelectField, IntegerField
import os
import configparser
from source import sequence_distribution, visualization
import json
import plotly
import plotly.graph_objects as go
from numpyencoder import NumpyEncoder
from plotly.offline import iplot
import shutil
from celery import Celery, Task
from celery import shared_task
from time import sleep
from celery.result import AsyncResult
import uuid
from pathlib import Path
from flask_mail import Mail, Message
import os


def celery_init_app(app: Flask) -> Celery:
    class FlaskTask(Task):
        def __call__(self, *args: object, **kwargs: object) -> object:
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config["CELERY"])
    celery_app.set_default()
    app.extensions["celery"] = celery_app
    return celery_app

app = Flask(__name__)
app.secret_key = 'tO$&!|0wkamvVia0?n$NqIRVWOG'
app.config['UPLOAD_FOLDER'] = 'static/sessions/'
app.config.from_mapping(
    CELERY=dict(
        broker_url="redis://localhost:6379/0",
        result_backend="redis://localhost:6379/0",
        task_ignore_result=True,
    ),
)
# Bootstrap-Flask requires this line
bootstrap = Bootstrap5(app)
# Flask-WTF requires this line
csrf = CSRFProtect(app)
celery_app = celery_init_app(app)


# Configuration for Flask-Mail
app.config['MAIL_SERVER'] = 'smtp.yandex.ru'
app.config['MAIL_PORT'] = 465
app.config['MAIL_USERNAME'] = os.environ.get('MAIL_USERNAME')
app.config['MAIL_PASSWORD'] = os.environ.get('MAIL_PASSWORD')
app.config['MAIL_USE_TLS'] = False
app.config['MAIL_USE_SSL'] = True


import secrets
foo = secrets.token_urlsafe(16)
app.secret_key = foo

# Initialize Flask-Mail
mail = Mail(app)

class SequenceItem(Form):
    type = StringField('Type')
    sequence = StringField('Sequence')

class InputForm(FlaskForm):
    session_name = StringField('Session Name', validators=[DataRequired()])
    items = FieldList(FormField(SequenceItem), min_entries=1, max_entries=10)
    limit = IntegerField('Limit', default=0)
    threshold = FloatField('Threshold', validators=[DataRequired()], default=0.9)
    smoothing = SelectField('Smoothing',
                choices=[('None','None'),('LOWESS','lowess'),('Whittaker Smoother', 'whittaker'),('savgol','savgol'),('confsmooth','confsmooth')])
    email = StringField('Email', validators=[DataRequired()])
    file = FileField('fastq_file', validators=[FileRequired()])
    submit = SubmitField('Submit')

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ['fastq']

def read_config(directory_path):
    config_file = os.path.join(directory_path, 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_file)
    parameters = dict(config['Parameters'])
    sequences = dict(config['Sequences'])
    return parameters, sequences

def process_directories(base_dir):
    all_sessions = []
    print(os.listdir(base_dir))
    for directory in os.listdir(base_dir):
        if os.path.isdir(os.path.join(base_dir, directory)):
            if os.path.isfile(os.path.join(base_dir, directory, 'config.ini')):
                session_id = directory
                parameters, sequences = read_config(os.path.join(base_dir, directory))
                imageURL = os.path.join(base_dir, directory, 'distribution.png') if \
                    os.path.isfile(os.path.join(base_dir, directory, 'distribution.png')) \
                    else None
                session_dict = {'sessionID': session_id,
                                'parameters': parameters,
                                'sequences': sequences,
                                'imageURL': imageURL}
                all_sessions.append(session_dict)
    return all_sessions

def gerenate_config(sequences,
                    parameters):
    config = configparser.ConfigParser()
    config['Sequences'] = {item['type']:item['sequence'] for item in sequences}
    config['Parameters'] = {'limit': parameters['limit'],
                            'threshold': parameters['threshold'],
                            'input_file': os.path.join(parameters['new_dir'], parameters['filename']),
                            'smoothing': parameters['smoothing']
                            }
    with open(os.path.join(parameters['new_dir'], 'config.ini'), 'w') as configfile:
        config.write(configfile)

@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm()
    if form.validate_on_submit():
        sequences = []
        for field in form.items.data:
            sequences.append({'type': field['type'],
                              'sequence': field['sequence']})
        limit = form.limit.data
        threshold = form.threshold.data
        smoothing = dict(form.smoothing.choices).get(form.smoothing.data)
        email = form.email.data

        file = form.file.data
        filename = secure_filename(file.filename)
        filename_noext = os.path.splitext(os.path.basename(filename))[0]
        # create a new directory named as a file name without an extension
        # new_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), app.config['UPLOAD_FOLDER'],
        #                        filename_noext)
        new_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), app.config['UPLOAD_FOLDER'],
                               str(form.session_name.data))
        try:
            os.mkdir(new_dir)
        except Exception as e:
            pass
        # save fastq file to the local server
        file.save(os.path.join(new_dir, filename))
        # session['session'] = new_dir
        session['session'] = form.session_name.data
        parameters = {'limit': limit,
                     'threshold': threshold,
                     'filename': filename,
                     'new_dir': new_dir,
                     'smoothing': smoothing,
                     'email': email,
                      'datetime': str(datetime.now())}

        gerenate_config(sequences,
                        parameters)

        input_data = {'sequences': {item['type']: item['sequence'] for item in sequences},
                      'parameters': parameters,
                      'session_name': form.session_name.data}
        session['input_data'] = input_data
        return redirect(url_for('results'))
    else:
        return render_template('index.html', form=form, page='index')


@app.route('/contacts')
def contacts():
    return render_template('contacts.html')

@app.route('/sessions')
def sessions():
    base_directory = app.config['UPLOAD_FOLDER']
    #print(base_directory)
    all_sessions_list = process_directories(base_directory)
    #
    # for session_info in all_sessions_list:
    #     print(session_info)

    return render_template('sessions.html', sessions=all_sessions_list, page='sessions')
#
@app.route('/experiment/<sessionID>')
def experiment(sessionID):
    base_directory = app.config['UPLOAD_FOLDER']
    directory_path = os.path.join(base_directory, sessionID)
    parameters, sequences = read_config(directory_path)
    parameters['filename'] = os.path.basename(parameters['input_file'])
    parameters['new_dir'] = directory_path
    data = {'sequences': sequences,
            'parameters': parameters,
            'session_name': sessionID}
    session['input_data'] = data
    return redirect(url_for('results'))

@app.route('/delete/<sessionID>')
def delete(sessionID):
    base_directory = app.config['UPLOAD_FOLDER']
    directory_path = os.path.join(base_directory, sessionID)
    try:
        shutil.rmtree(directory_path)
        print(f"Directory '{directory_path}' and all its contents have been successfully removed.")
    except Exception as e:
        print(f"An error occurred while removing '{directory_path}': {e}")
    session['input_data'] = None
    return redirect(url_for('sessions'))

def generate_session_id():
    # Generate a random UUID (version 4)
    session_id = uuid.uuid4()
    return str(session_id)

@shared_task(ignore_result=False)
def data_processing(data):
    output_data = sequence_distribution.main(data['parameters']['new_dir'])
    with open(os.path.join(data['parameters']['new_dir'], 'sequences.json'), 'w') as f:
        json.dump(output_data, f, default=str)

    fig1 = visualization.plot_distribution_proportions(output_data['sequences'], data['parameters']['smoothing'])
    distrJSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)

    with open(os.path.join(data['parameters']['new_dir'], 'distribution.png'), "wb") as distribution_file:
        fig1.write_image(distribution_file)

    return {'session_id': Path(data['parameters']['new_dir']).name,
            'email': data['parameters']['email']}

@app.route('/results')
def results():
    if 'input_data' in session and session['input_data'] is not None:
        data = session['input_data']
        if os.path.exists(os.path.join(data['parameters']['new_dir'], 'sequences.json')):
            with open(os.path.join(data['parameters']['new_dir'], 'sequences.json'), 'r') as f1:
                output_data = json.load(f1)

                fig1 = visualization.plot_distribution_proportions(output_data['sequences'], data['parameters']['smoothing'])
                distrJSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)
                sequences = [{'type': seq['type'],
                              'sequence': seq['sequence'],
                              'peaks': seq['peaks'] if 'peaks' in seq else [],
                              'noise_level': seq['noise_level'] if 'noise_level' in seq else 0,
                              'total_reads': seq['total_reads'] if 'total_reads' in seq else 0,
                              'total_proportion': seq['total_proportion'] if 'total_proportion' in seq else 0,
                              'value_counts': seq['value_counts'] if 'value_counts' in seq else []
                              } for seq in output_data['sequences']]
                fastq_parameters = {'n_records': output_data['parameters']['n_records'],
                                    'avg_noise_level': output_data['parameters']['avg_noise_level']}

                return render_template('results.html',
                                       plots={'hist1': distrJSON},
                                       data=data,
                                       sequences=sequences,
                                       fastq_parameters=fastq_parameters,
                                       page='results')
        else:
            result = data_processing.delay(data)
            return render_template('async_result.html',
                                   result_id=result.id,
                                   parameters = data['parameters'],
                                   page='results')
    else:
        return render_template('no_results.html', page='results')

def send_email(email, sessionID, path):
    recipient = email
    subject = f'NanoporeInspect: session {sessionID} results are ready'
    message_body = f'The results are available in the web application by the link: {path}'
    msg = Message(
        subject=subject,
        sender=app.config['MAIL_USERNAME'],
        recipients=[recipient]
    )
    msg.body = message_body  # Plain text email body

    try:
        mail.send(msg)
        return f"Email sent to {recipient}!"
    except Exception as e:
        return f"Failed to send email. Error: {str(e)}"

@app.route('/send_email')
def send_email_test():
    recipient = 'magsend@gmail.com'
    subject = 'TEST from NanoporeInspect'
    message_body = 'Test from NanoporeInspect'

    msg = Message(
        subject=subject,
        sender=app.config['MAIL_USERNAME'],
        recipients=[recipient]
    )
    msg.body = message_body  # Plain text email body

    try:
        mail.send(msg)
        return f"Email sent to {recipient}!"
    except Exception as e:
        return f"Failed to send email. Error: {str(e)}"

@app.route("/result/<id>", methods=['GET','POST'])
def task_result(id: str) -> dict[str, object]:
    result = AsyncResult(id)
    if result.ready():
        session_id = result.result['session_id']
        email = result.result['email']
        path = request.url_root + 'experiment/' + session_id
        send_email(email, session_id, path)
        return redirect(url_for('experiment', sessionID=session_id))
    else:
        return render_template('in_progress.html', result_id = id, page='results')

    # return {
    #     "ready": result.ready(),
    #     "successful": result.successful(),
    #     "value": result.result if result.ready() else None,
    # }

if __name__ == '__main__':
    app.run(debug=True)