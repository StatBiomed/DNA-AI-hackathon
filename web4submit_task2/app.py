from flask import Flask, render_template, request, jsonify
import os
import pandas as pd
import subprocess

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = "results"
# RESULTS_FILE = 'results.csv'

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/results/<task_id>', methods=['GET'])
def get_results(task_id):
    # Load results from the CSV file
    results_file = f"{RESULTS_FOLDER}/results_{task_id}.csv"
    if os.path.exists(results_file):
        try:
            results_df = pd.read_csv(results_file).fillna("")
            results_json = results_df.to_dict(orient='records')
            return jsonify(results_json)
        except Exception as e:
            return jsonify({'error': str(e)}), 500
    else:
        return jsonify([])

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files or 'modelName' not in request.form or 'taskId' not in request.form:
        return jsonify({'error': 'File or Model Name or Task missing'}), 400

    file = request.files['file']
    model_name = request.form['modelName']
    taskId = request.form['taskId']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
    file.save(filepath)

    try:
        # Execute Python script on the uploaded file
        #print(filepath)
        result = subprocess.check_output(['python', 'process_file.py', filepath, taskId], text=True)

        #print(result)
        # Create a DataFrame for the result
        result_data = result.strip().split('\n')
        #print(result_data)
        result_df = pd.DataFrame([row.split(',') for row in result_data[1:]], columns=result_data[0].split(','))

        # Add the model name as a new column
        result_df['Model Name'] = model_name

        # Append or create the results CSV
        results_file = f"{RESULTS_FOLDER}/results_{taskId}.csv"
        if os.path.exists(results_file):
            existing_df = pd.read_csv(results_file)
            updated_df = pd.concat([existing_df, result_df])
            updated_df.to_csv(results_file, index=False)
        else:
            result_df.to_csv(results_file, index=False)

        # Return processed data for immediate display
        return jsonify({'data': result})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':

    app.run(host='10.64.155.14', port=5013, debug=True)
    ## test on localhost
    # app.run(host='localhost', port=5011, debug=True)
