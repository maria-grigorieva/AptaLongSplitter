import pandas as pd
from flask import Flask, render_template, request, url_for, redirect, session
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


app = Flask(__name__)
app.secret_key = 'tO$&!|0wkamvVia0?n$NqIRVWOG'
app.config['UPLOAD_FOLDER'] = 'static/sessions/'
# Bootstrap-Flask requires this line
bootstrap = Bootstrap5(app)
# Flask-WTF requires this line
csrf = CSRFProtect(app)

import secrets
foo = secrets.token_urlsafe(16)
app.secret_key = foo

class SequenceItem(Form):
    type = StringField('Type')
    sequence = StringField('Sequence')

class InputForm(FlaskForm):
    session_name = StringField('Session Name', validators=[DataRequired()])
    items = FieldList(FormField(SequenceItem), min_entries=1, max_entries=10)
    limit = IntegerField('Limit', default=500)
    threshold = FloatField('Threshold', validators=[DataRequired()], default=0.75)
    smoothing = SelectField('Smoothing',
                choices=[('None','None'),('LOWESS','lowess'),('Whittaker Smoother', 'whittaker'),('savgol','savgol'),('confsmooth','confsmooth')])
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
                     'smoothing': smoothing}

        gerenate_config(sequences,
                        parameters)

        input_data = {'sequences': {item['type']: item['sequence'] for item in sequences},
                      'parameters': parameters,
                      'session_name': form.session_name.data}
        session['input_data'] = input_data
        return redirect(url_for('results'))
    else:
        return render_template('index.html', form=form)


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

    return render_template('sessions.html', sessions=all_sessions_list)
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

# def plotlyfromjson(json_plot):
#     fig = go.Figure(data=json_plot['data'], layout=json_plot['layout'])
#     return fig


@app.route('/results')
def results():
    if 'input_data' in session:
        data = session['input_data']
        if os.path.exists(os.path.join(data['parameters']['new_dir'], 'sequences.json')):
            with open(os.path.join(data['parameters']['new_dir'], 'sequences.json'), 'r') as f1:
                output_data = json.load(f1)

                fig1 = visualization.plot_distribution_proportions(output_data['sequences'], data['parameters']['smoothing'])
                distrJSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)
                #
                # peaks_table = visualization.make_peaks_subplots(sequences)
                # peaksJSON = json.dumps(peaks_table, cls=plotly.utils.PlotlyJSONEncoder)
                sequences = [{'type': seq['type'],
                              'sequence': seq['sequence'],
                              'peaks': seq['peaks'],
                              'noise_level': seq['noise_level'],
                              'total_reads': seq['total_reads'],
                              'total_proportion': seq['total_proportion']
                              } for seq in output_data['sequences']]
                fastq_parameters = {'n_records': output_data['parameters']['n_records'],
                                    'avg_noise_level': output_data['parameters']['avg_noise_level']}

                return render_template('results.html',
                                       plots={'hist1': distrJSON},
                                       data=data,
                                       sequences=sequences,
                                       fastq_parameters=fastq_parameters)
        else:
            output_data = sequence_distribution.main(data['parameters']['new_dir'])
            with open(os.path.join(data['parameters']['new_dir'],'sequences.json'), 'w') as f:
                json.dump(output_data, f, default=str)

            fig1 = visualization.plot_distribution_proportions(output_data['sequences'], data['parameters']['smoothing'])
            distrJSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)

            with open(os.path.join(data['parameters']['new_dir'],'distribution.png'), "wb") as distribution_file:
                fig1.write_image(distribution_file)

            plots = {'hist1': distrJSON}
            sequences = [{'type': seq['type'],
                          'sequence': seq['sequence'],
                          'peaks': seq['peaks'],
                          'noise_level': seq['noise_level'],
                          'total_reads': seq['total_reads'],
                          'total_proportion': seq['total_proportion']
                          } for seq in output_data['sequences']]
            fastq_parameters = {'n_records': output_data['parameters']['n_records'],
                                'avg_noise_level': output_data['parameters']['avg_noise_level']}
            return render_template('results.html',
                                   plots=plots,
                                   data=data,
                                   sequences=sequences,
                                   fastq_parameters=fastq_parameters)
    else:
        return render_template('no_results.html')

if __name__ == '__main__':
    app.run(debug=True)