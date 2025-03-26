from flask import Flask, render_template, request, jsonify
import os
import pandas as pd
import subprocess

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
RESULTS_FILE = 'results.csv'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/results', methods=['GET'])
def get_results():
    # Load results from the CSV file
    if os.path.exists(RESULTS_FILE):
        try:
            results_df = pd.read_csv(RESULTS_FILE)
            results_json = results_df.to_dict(orient='records')
            return jsonify(results_json)
        except Exception as e:
            return jsonify({'error': str(e)}), 500
    else:
        return jsonify([])

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files or 'modelName' not in request.form:
        return jsonify({'error': 'File or Model Name missing'}), 400

    file = request.files['file']
    model_name = request.form['modelName']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
    file.save(filepath)

    try:
        # Execute Python script on the uploaded file
        result = subprocess.check_output(['python', 'process_file.py', filepath], text=True)

        # Create a DataFrame for the result
        result_data = result.strip().split('\n')
        result_df = pd.DataFrame([row.split(',') for row in result_data[1:]], columns=result_data[0].split(','))

        # Add the model name as a new column
        result_df['Model Name'] = model_name

        # Append or create the results CSV
        if os.path.exists(RESULTS_FILE):
            existing_df = pd.read_csv(RESULTS_FILE)
            updated_df = pd.concat([existing_df, result_df])
            updated_df.to_csv(RESULTS_FILE, index=False)
        else:
            result_df.to_csv(RESULTS_FILE, index=False)

        # Return processed data for immediate display
        return jsonify({'data': result})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='10.64.155.14', port=5011, debug=True)
