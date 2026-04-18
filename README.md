# SWGDRUG Library Parser

## Description
A Python utility designed to parse the standard SWGDRUG `.txt` library, export all compound metadata into a consolidated CSV, and generate styled, dark-themed mass spectrum plots for every compound in the library. This script is ideal for translating rigid text-based spectral libraries into highly visual, accessible formats for forensic reference and presentations.

## Features
* **Automated Parsing:** Reads the SWGDRUG text file and extracts key metadata for each compound, including Name, Formula, MW, ExactMass, CASNO, ID, and Comment.
* **CSV Export:** Compiles all parsed compound metadata, along with the calculated total number of peaks (`Num_Peaks`), into a single `SWGDRUG_compounds.csv` file.
* **High-Fidelity Visualizations:** Generates an individual `.png` mass spectrum plot for each compound. The plots feature a clean dark theme (background `#0d1117`), subtle gridlines, and dynamic bar coloring where taller bars receive lighter blue tones.
* **Smart Annotation:** Automatically identifies and labels the top 10 peaks with their *m/z* values based on their relative abundance. Relative abundance is normalized to a 0-999 scale for consistency.
* **Sanitized Outputs:** Automatically cleans filenames by removing invalid characters (`\ / * ? : " < > |`) to prevent operating system errors when saving the image files.

## Requirements
This script relies on external data manipulation and plotting libraries. You will need to install them using pip:
`pip install matplotlib pandas`

## Configuration & Usage
1. Open the script in your preferred editor.
2. Update the configuration variables at the top of the file to point to your specific directories:
   * `INPUT_FILE`: The full file path to your downloaded SWGDRUG `.txt` file (e.g., `r"C:\Users\txbar\Downloads\SWGDRUG 314.txt"`).
   * `OUTPUT_DIR`: The folder where you want the CSV and all the PNG spectrum plots to be saved.
3. Run the script via your terminal or IDE:
   `python swgdrug_parser.py`
4. The script will print its progress to the console, notifying you as it parses the compounds, writes the CSV, and generates each plot.
